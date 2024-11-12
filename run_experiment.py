''' Volume-Controlled Cavitation Data Collection System (VCC-DCS)
Authors: Joseph Beckett, Brandon Brozich, and Esben Nielsen 
Description: Script for running needle-induced cavitation experiments in a custom experimental set-up
             developed at the University of Michigan under Prof. Jon Estrada starting in 2024. The script
             manages test definition, data acquisition, and data export to .csv files for additional 
             post-processing and inverse calibration of hyperelastic constituitive parameters in MATLAB.'''


#### INPUTS SECTION ####
## Troubleshooting Commands: Set to <True> to print live volume, pressure, and force readings to the terminal without pumping fluid.
# Note: To exit trouble shotting mode, click on the terminal where live readings are being printed and hit Ctrl+C.
MODE_LIVE_READ = False

## System Defintions: Define parameters that only change with alterations to the physical equiptment.
DIAMETER_SYRINGE = 1.03     # Inner diameter of syringe (Hamilton Model 1705) in mm.
CAPACITY_SYRINGE = 50       # Syringe (Hamilton Model 1705) volumetric capacity in microliters (µL).
GAIN_FORCE_AMP = 1015.66    # Should always be 1015.66 unless DIP switches in the load cell Futek amplifier are changed from [0 0 0 0 1 1 1 0]
GAIN_PRESSURE_AMP = 106.26  # Should always be 106.26 unless DIP switches in the load cell Futek amplifier are changed from [0 0 0 0 0 0 0 1]
VOLTAGE_RANGE_PRESSURE = [-10.0,10.0] # Voltage range (V) of pressure readings entering DAQ. Options limited to ±0.2 V, ±1 V, ±5 V, ±10 V.
VOLTAGE_RANGE_FORCE = [-10.0,10.0]   
LOAD_RANGE_PRESSURE = [-68.9,413.6] # (kPa) ~ (-10 to 60 psi) Pressre extrema for Pendotech PRESS-S-000.
LOAD_RANGE_FORCE = [-222.4,222.4] # (N) ~ (±50 lb) Force extrema for Futek LCM300.

## Experiment Settings
VOLUME_INJECTION = 881.8    # (nL) Amount of volume to be injected.
FLOW_RATE = 200             # (nL/s) Constant flow rate.

# Data Aquisition Settings
# Note: If analog input (i.e., force and pressure) readings are collected in excess of 1 kHz, ensure the 'Bandwidth' setting in the amplifier is set accordingly.
SAMPLE_RATE = 200           # (Hz) Sampling rate of all channels.                     
WAIT_TIME_PREINJECTION = 1  # (s) Duration that force and pressure readings are collected prior to starting the injection.
WAIT_TIME_POSTINJECTION = 1 # (s) Duration that force and pressure readings are collected after the injection is completed.
TIMED_DAQ_END = True     # If True, data aquisition ends after WAIT_TIME_POSTINJECTION, otherwise type a character in terminal and hit ENTER to manually 
                            # stop at any time after start of injection
from datetime import datetime
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # Format: YYYYMMDD_HHMMSS
FOLDERPATH = 'C:\\Users\\jbecktt\\Desktop\\VCC_System\\VCC_test_'+timestamp+'\\' # This folder will be made if it does not already exist.
FILEPATH_CSV = FOLDERPATH + 'VCC_'+timestamp+'.csv'

#### CLASS DEFINITIONS ####
import serial, time, math, csv, nidaqmx, serial.tools.list_ports, os, sys, numpy as np, _thread
from nidaqmx.constants import AcquisitionType, TerminalConfiguration, EncoderType, AngleUnits, Edge, Level, FrequencyUnits, LoggingOperation, LoggingMode, READ_ALL_AVAILABLE
from nptdms import TdmsFile

class TDMSHandler:
    def __init__(self, folderpath,quiet=False):
        self.folderpath = folderpath
        self.quiet = quiet
        self.chkdir()

    def chkdir(self):
        if not os.path.exists(self.folderpath):
            try:
                os.makedirs(self.folderpath)
                if not self.quiet:
                    print(f"Directory ({self.folderpath}) did not exist. New directory has been made.")
            except:
                raise Exception("Error making directory {self.folderpath}.")
            
    def chkfile(filepath):
        if os.path.exists(filepath):
            raise Exception("File ({filepath}) already exists. Test stopped to prevent overwriting data.")

    def read_tdms_data(meas_obj,raw=False):
        with TdmsFile.read(meas_obj.filepath) as tdms_file:
            if len(tdms_file.groups()) != 1: #or len(tdms_file.groups()[0].channels()) != 1:
                raise Exception("TDMS file must contain exactly one group with exactly one channel. Code must be altered to accept more general data.")
            if isinstance(meas_obj, AnalogInput):
                data = {}
                for i in range(len(meas_obj.ai_chan)):
                    voltage = tdms_file["analog_in"][meas_obj.ai_chan[i]].data
                    if raw == True:
                        data[i] = voltage #(V)
                    else:
                        if meas_obj.data_type[i] == "pressure":
                            data[i] = meas_obj._voltageV_to_pressurekPa(voltage) # (kPa)
                        elif meas_obj.data_type[i] == "force":
                            data[i] = meas_obj._voltageV_to_forceN(voltage) # (N)
                        else: raise Exception("Data type not supported yet for .tdms file reading.")
                return data

            # if isinstance(meas_obj, PressureSensor):
            #     voltage = tdms_file["pressure"][meas_obj.ai_chan].data # V
            #     pressure = meas_obj._voltageV_to_pressurekPa(voltage) # kPa
            #     return pressure if not raw else voltage
            # elif isinstance(meas_obj, LoadCell):
            #     voltage = tdms_file["force"][meas_obj.ai_chan].data # V
            #     force = meas_obj._voltageV_to_forceN(voltage) # N
            #     return force if not raw else voltage
            elif isinstance(meas_obj, OpticalEncoder):
                cnt = tdms_file["angular_disp"][meas_obj.ci_chan].raw_data # count
                position_incr = tdms_file["angular_disp"][meas_obj.ci_chan].properties["NI_Scale[1]_PositionEncoder_Position_Increment"] # angular_disp/count
                angle = cnt*position_incr # degree
                volume = meas_obj._angle_deg_to_vol_nL(angle) # nL
                return volume if not raw else cnt
            else:
                raise Exception("Measurement object type not supported yet for .tdms file reading.")

class NIDAQTask:
    def __init__(self):
        self.task = nidaqmx.Task()

    def start(self):
        return self.task.start()
    
    def stop(self):
        return self.task.stop()

    def close(self):
        self.task.close()

    def commit_task(self):
        self.task.control(nidaqmx.constants.TaskMode.TASK_COMMIT) # Pre-commit some resources and settings prior to starting task to reduce lag of task.start() command

    def read(self):
        return self.task.read()   
    
    def read_all(self):
        return self.task.read(number_of_samples_per_channel=READ_ALL_AVAILABLE) # If sample mode is CONTINUOUS, all data points currently buffer are read.
                                                                                # If sample mode is FINITE, the code waits until the task is finished and then reads the buffer.

class OpticalEncoder(NIDAQTask):
    def __init__(self, ci_chan, decoding_type, diameter_syringe, filepath, pulse_gen=None, sample_rate=None, n_samples=None, sample_mode=AcquisitionType.CONTINUOUS):
        super().__init__()
        units = AngleUnits.DEGREES
        pulses_per_rev = 100 # E4T encoder spec
        self.task.ci_channels.add_ci_ang_encoder_chan(counter=ci_chan,
                                                      decoding_type=decoding_type,
                                                      pulses_per_rev=pulses_per_rev,
                                                      units=units,
                                                      initial_angle=0.0,
                                                      zidx_enable=False # E4T encoder does not have z-index channel
                                                      )
        if isinstance(pulse_gen,PulseGenerator): # Set time source to pulse generator (if provided)
            self.task.timing.cfg_samp_clk_timing(rate=pulse_gen.sample_rate,
                                                 source=pulse_gen.co_pulse_term,
                                                 active_edge=Edge.RISING,
                                                 sample_mode=pulse_gen.sample_mode,
                                                 samps_per_chan=pulse_gen.samps_per_chan)
        elif sample_rate is not None:
            if n_samples is not None:
                self.task.timing.cfg_samp_clk_timing(rate=sample_rate,
                                        sample_mode=sample_mode,
                                        samps_per_chan=n_samples) # Including n_samples is necessary for AcquisitionType.FINITE, but it is helpful for estimating memory usage in AcquisitionType.CONTINUOUS mode
            else:
                self.task.timing.cfg_samp_clk_timing(sample_rate=sample_rate,
                                        sample_mode=sample_mode)
        if isinstance(pulse_gen,PulseGenerator) or sample_rate is not None:
            TDMSHandler.chkfile(filepath)
            self.task.in_stream.configure_logging(filepath,
                                                  LoggingMode.LOG_AND_READ,    
                                                  group_name='angular_disp',
                                                  operation=LoggingOperation.CREATE) # Create a new TDMS file. If the file already exists, NI-DAQmx returns an error.
        self.ci_chan = ci_chan
        self.diameter_syringe = diameter_syringe
        self.filepath = filepath

    def _angle_deg_to_vol_nL(self,angle_deg):
        gear_ratio = 1/2.4
        factor = gear_ratio*(1/24)*(25400/1)*(1/360)*10**-3 # mm/°of spindle rotation
        A_syringe = math.pi*(self.diameter_syringe/2)**2 # mm^2
        vol = angle_deg*factor*A_syringe # mm^3 or μL
        return vol*10**3 # nL

    def read_volume(self):
        return self._angle_deg_to_vol_nL(self.read())
    
    def print_volume(self):
        print("%.3f nL" % self.read_volume())

class AnalogInput(NIDAQTask):
    def __init__(self, ai_chan, data_type, gain, ai_range_V, ai_range_conv, filepath, pulse_gen=None, sample_rate=None, n_samples=None, sample_mode=AcquisitionType.CONTINUOUS):
        super().__init__() # Ensures that NIDAQTask Class is also initialized
        for i in range(len(ai_chan)):
            self.task.ai_channels.add_ai_voltage_chan(physical_channel=ai_chan[i], 
                                                      terminal_config=TerminalConfiguration.DIFF, # reflects currently wiring (decision informed by NI best practices)
                                                      min_val = ai_range_V[i][0], max_val = ai_range_V[i][1])
        if isinstance(pulse_gen,PulseGenerator): # Set time source to pulse generator (if provided)
            self.task.timing.cfg_samp_clk_timing(rate=pulse_gen.sample_rate, 
                                                    source=pulse_gen.co_pulse_term, 
                                                    active_edge=Edge.RISING, 
                                                    sample_mode=pulse_gen.sample_mode)
        elif sample_rate is not None:
            if n_samples is not None:
                self.task.timing.cfg_samp_clk_timing(rate=sample_rate,
                                        sample_mode=sample_mode,
                                        samps_per_chan=n_samples) # Including n_samples is necessary for AcquisitionType.FINITE, but it is helpful for estimating memory usage in AcquisitionType.CONTINUOUS mode
            else:
                self.task.timing.cfg_samp_clk_timing(sample_rate=sample_rate,
                                        sample_mode=sample_mode)
        if isinstance(pulse_gen,PulseGenerator) or sample_rate is not None:
            TDMSHandler.chkfile(filepath)
            self.task.in_stream.configure_logging(filepath,
                                                  LoggingMode.LOG_AND_READ,
                                                  group_name='analog_in',
                                                  operation=LoggingOperation.CREATE) # Create a new TDMS file. If the file already exists, NI-DAQmx returns an error.
        self.ai_chan = ai_chan
        self.data_type = data_type
        self.gain = gain
        self.ai_range_V = ai_range_V
        self.ai_range_conv = ai_range_conv
        self.filepath = filepath
        self._get_overvoltage_limits()

    def _voltageV_to_pressurekPa(self,voltage):
        mV_psi_Vexic = 0.2584 # mV/psi
        Vexcit = 10 # excitation voltage of pressure sensor (V)
        kPa_psi = 6.8947572932 # kPa/psi
        mV_kPa = mV_psi_Vexic*Vexcit/kPa_psi
        pressure = voltage*10**3/(self.gain[self.data_type.index("pressure")]*mV_kPa) # kPa
        return pressure
    
    def read_pressure(self,raw = False):
        voltage = self.read()
        return self._voltageV_to_pressurekPa(voltage[self.data_type.index("pressure")]) if not raw else voltage[self.data_type.index("pressure")] # only print the most recent value

    def print_pressure(self):
        print("%.3f kPa" % self.read_pressure())

    def _voltageV_to_forceN(self,voltage):
        N_lb = 4.44822 # N/lb conversion
        FSR = 50 # full-scale-range on loadcell (lb)
        Vexcit = 10 # excitation voltage of loadcell (V)
        RO = 1 # rated output of load cell (mV/V)
        sensitivity = FSR/(Vexcit*RO*10**-3*self.gain[self.data_type.index("force")]) # lb/V
        force = sensitivity*voltage*N_lb # N
        return force
    
    def read_force(self, raw = False):
        voltage = self.read()
        return self._voltageV_to_forceN(voltage[self.data_type.index("force")]) if not raw else voltage[self.data_type.index("force")] # only print the most recent value

    def print_force(self):
        print("%.3f N" % self.read_force())
    
    def _get_overvoltage_limits(self):
        self.overvoltage_range = {}; self.overload = {}
        for i in range(len(self.ai_chan)):
            if self.data_type[i] == "force":
                Vlim_low = [self.ai_range_V[i][0],self.ai_range_conv[i][0]/self._voltageV_to_forceN(1)]
                Vlim_high = [self.ai_range_V[i][1],self.ai_range_conv[i][1]/self._voltageV_to_forceN(1)]
            elif self.data_type[i] == "pressure":
                Vlim_low = [self.ai_range_V[i][0],self.ai_range_conv[i][0]/self._voltageV_to_pressurekPa(1)]
                Vlim_high = [self.ai_range_V[i][1],self.ai_range_conv[i][1]/self._voltageV_to_pressurekPa(1)]
            else: raise Exception("Data type not supported yet for overvoltage protection.")
            self.overvoltage_range[i] = [max(Vlim_low), min(Vlim_high)]
            self.overload[i] = [Vlim_low.index(self.overvoltage_range[i][0]), Vlim_high.index(self.overvoltage_range[i][1])]

    def safe_voltage_check(self,data):
        data = np.array(data)
        for i in range(len(self.ai_chan)):
            overvoltage = np.any(data[i] < self.overvoltage_range[i][0]) or np.any(data[i] > self.overvoltage_range[i][1])
            if overvoltage == True:
                return False
            else : return True
    
class PulseGenerator(NIDAQTask):
    def __init__(self,co_pulse_term,co_chan,sample_rate,n_samples=None,idle_state=Level.LOW,initial_delay=0.0,duty_cycle=0.5):
        super().__init__()
        sample_mode = AcquisitionType.CONTINUOUS
        co_channel = self.task.co_channels.add_co_pulse_chan_freq(co_chan, idle_state=idle_state, 
                                                             initial_delay=initial_delay, 
                                                             freq=sample_rate, 
                                                             duty_cycle=duty_cycle, 
                                                             units=FrequencyUnits.HZ)
        co_channel.co_pulse_term = co_pulse_term # Output pulse train on a specified terminal (PFI line)
        self.task.timing.cfg_implicit_timing(sample_mode=sample_mode)
        self.co_pulse_term = co_pulse_term
        self.sample_rate = sample_rate
        self.samps_per_chan = n_samples
        self.sample_mode = sample_mode

class SyringePump:
    def query(self, command):
        self.serial_port.write((command + '\r').encode())
        response = bytearray()
        start_time = time.time()
        while time.time() - start_time < 5: # Time out after 5 seconds
            bytes_waiting = self.serial_port.in_waiting
            if bytes_waiting:
                chunk = self.serial_port.read(bytes_waiting)
                response.extend(chunk)
                if b'\x11' in response:
                    index = response.index(b'\x11')
                    return response[:index].decode()
        raise Exception(f"Unable to query syringe pump. Timed out with response: {response.decode()}")

    
    def __init__(self, port):
        self.serial_port = serial.Serial(port, 115200, timeout=1) # Opens serial connection with syringe pump (~3.5 ms response time)
        self.query('poll on') # 'poll on' must be the first query sent to the syringe pump.

    def is_moving(self):
        stat = self.query('@status')
        if stat[-1] == '>':
            return True
        elif stat[-1] == '<':
            return True
        else: 
            return False

    def get_infused_vol(self):
        '''Returns infused volume in nanoliters (nL).'''
        status_response = self.query('@status')
        status_fields = status_response.split()
        infused_volume = int(status_fields[2])
        return infused_volume*1e-6
    
    def set_flow_rate(self, flow_rate):
        '''Sets flow rate of pump to flow_rate (nL/s).'''
        if self.query('status').split()[3][0].capitalize() == 'I':
            direction = 'i'
        elif self.query('status').split()[3][0].capitalize() == 'W':
            direction = 'w'
        else:
            raise Exception('Flow rate could not be set:\n Pump is not in infuse or withdraw mode.')
        start = time.perf_counter()
        while not math.isclose(float(self.query('@' + direction + 'rate').split()[0]), flow_rate):
            set_rate_msg = self.query('@' + direction + 'rate ' + str(flow_rate) + ' nL/s')
            if 'error' in set_rate_msg:
                raise Exception('Flow rate could not be set:\n' + set_rate_msg)
            if time.perf_counter() - start > 5: # 5 second timeout
                raise Exception('Flow rate could not be set: Timed out.')
            time.sleep(0.01)

    def infuse_const_rate(self, target_vol, flow_rate):
        '''Drives syringe forward up to target_vol (nL) at constant rate of flow_rate (nL/s).\n
        WARNING: Starting the pump takes a while (sometimes up to a second). Avoid calling this function often.\n
        To change flow rate mid-withdraw, call set_flow_rate().'''
        load_msg = self.query('load qs i')
        if 'error' in load_msg:
            raise Exception('Infuse program not loaded:\n' + load_msg)
        start_load = time.perf_counter()
        while self.query('@status').split()[3][0].capitalize() != 'I':
            if time.perf_counter() - start_load > 5: # 5 second timeout
                raise Exception('Failed to load quickstart program: Timed out.')
            time.sleep(0.01)
        set_vol_msg = self.query('tvolume ' + str(target_vol) + ' nL')
        if 'error' in set_vol_msg:
            raise Exception('Target volume could not be set:\n' + set_vol_msg)
        self.set_flow_rate(flow_rate)
        start_run = time.perf_counter()
        while time.perf_counter() - start_run < 5: # 5 second timeout # JGB temporary comment on 8/8/2024
            self.query('@run')
            time.sleep(0.01)
            if self.query('@status').split()[3][0] == 'I':
                return
        raise Exception('Failed to start pump.')

    def withdraw_const_rate(self, target_vol, flow_rate):
        '''Drives syringe in reverse up to target_vol (nL) at constant rate of flow_rate (nL/s).\n
        WARNING: Starting the pump takes a while (sometimes up to a second). Avoid calling this function often.\n
        To change flow rate mid-withdraw, call set_flow_rate().'''
        load_msg = self.query('load qs w')
        if 'error' in load_msg:
            raise Exception('Withdraw program not loaded:\n' + load_msg)
        start_load = time.perf_counter()
        while self.query('status').split()[3][0].capitalize() != 'W':
            if time.perf_counter() - start_load > 5: # 5 second timeout
                raise Exception('Failed to load quickstart program: Timed out.')
            time.sleep(0.01)
        set_vol_msg = self.query('tvolume ' + str(target_vol) + ' nL')
        if 'error' in set_vol_msg:
            raise Exception('Target volume could not be set:\n' + set_vol_msg)
        self.set_flow_rate(flow_rate)
        start_run = time.perf_counter()
        while time.perf_counter() - start_run < 5: # 5 second timeout
            self.query('@run')
            time.sleep(0.01)
            if self.query('@status').split()[3][0] == 'W':
                return
        raise Exception('Failed to start pump.')
    
    def stop(self):
        self.query('stop')
        status = self.query('@status')[-1]
        start = time.perf_counter()
        while status != '*' and status != ':':
            time.sleep(0.01)
            if time.perf_counter() - start > 5: # 5 second timeout
                raise Exception('Unable to stop syringe pump: Timed out.')
            status = self.query('@status')[-1]

    def set_syringe_dims(self, diameter, volume):
        '''Sets syringe dimensions: diameter (mm), volume (uL).\n
        Call this function before setting rates and target volume.'''
        start_time_diameter = time.perf_counter()
        while not math.isclose(float(self.query('diameter').split()[0]), diameter):
            self.query('diameter ' + str(diameter))
            time.sleep(0.01)
            if time.perf_counter() - start_time_diameter > 5: # 5 second timeout
                raise Exception('Unable to set syringe diameter. Timed out')

        start_time_volume = time.perf_counter()
        self.query('svolume ' + '0' + ' ul') # Setting units to ul
        while not math.isclose(float(self.query('svolume').split()[0]), volume):
            self.query('svolume ' + str(volume) + ' ul')
            time.sleep(0.01)
            if time.perf_counter() - start_time_volume > 5: # 5 second timeout
                raise Exception('Unable to set syringe volume. Timed out')
            
    def __del__(self):
        if self.serial_port.is_open: self.serial_port.close()

#### TEST SET-UP ####
if not MODE_LIVE_READ: TDMSHandler(FOLDERPATH)
global pump, analog_in, TEST_ACTIVE
## Initializing Syringe Pump and Sensors
pump = SyringePump('COM4')
pump.set_syringe_dims(DIAMETER_SYRINGE, CAPACITY_SYRINGE)
duration_injection = VOLUME_INJECTION/FLOW_RATE
time_DAQ_record = WAIT_TIME_PREINJECTION + duration_injection + WAIT_TIME_POSTINJECTION
n_samples = math.ceil(time_DAQ_record*SAMPLE_RATE)
pulse_gen = PulseGenerator('/Dev2/PFI12','/Dev2/ctr1',SAMPLE_RATE,n_samples) if not MODE_LIVE_READ else None
encoder_type = EncoderType.X_4 # X_2 is generally has lower resolution than X_4, but is less sensitive to vibrations
optical_encoder = OpticalEncoder('Dev2/ctr0', encoder_type, DIAMETER_SYRINGE, FOLDERPATH+'angular_disp.tdms', pulse_gen)
analog_in = AnalogInput(['Dev2/ai0','Dev2/ai1'], ["pressure","force"], [GAIN_PRESSURE_AMP, GAIN_FORCE_AMP], [VOLTAGE_RANGE_PRESSURE,VOLTAGE_RANGE_FORCE], [LOAD_RANGE_PRESSURE,LOAD_RANGE_FORCE], FOLDERPATH+'analog_in.tdms', pulse_gen)

# Start system protection thread to protect from overvoltage of DAQ or overload of sensor(s)
def overvoltage_protect_thread():
    global pump, analog_in, TEST_ACTIVE
    while TEST_ACTIVE:
        ai_data = analog_in.read() if MODE_LIVE_READ else analog_in.read() # the latter was changed from read_all() to read() during troubleshooting!!!
        if MODE_LIVE_READ == True: 
            time.sleep(0.001)
        safe_voltage = analog_in.safe_voltage_check(ai_data)
        if not safe_voltage:
            TEST_ACTIVE = False
            pump.stop()
            print("NIC Update: Pump stopped due to overvoltage.")

def print_live_sensor_readout(t_wait=0.1):
    try: 
        while MODE_LIVE_READ:
            vol_str = "Volume: %.3f nL  " % optical_encoder.read_volume()
            pressure_V = analog_in.read_pressure(raw=True)
            pressure_str = "Pressure: %.3f kPa (%.5f V)  " % (analog_in._voltageV_to_pressurekPa(pressure_V), pressure_V)
            force_V = analog_in.read_force(raw=True)
            force_str = "Force: %.3f N (%.5f V)  " % (analog_in._voltageV_to_forceN(force_V), force_V)
            print(f"{vol_str}\t{pressure_str}\t{force_str}", end="\r")
            time.sleep(t_wait) 
    except nidaqmx.DaqError as e:
        print('An error occurred during live reading mode: ', e)
    except KeyboardInterrupt:
        print("\nLive reading mode stopped by user")
    finally: 
        optical_encoder.stop(); optical_encoder.close()
        analog_in.stop(); analog_in.close()
    sys.exit() # Exit script after live reading mode is stopped

#### RUN VCC TEST ####
## Start Data Acquisition and Injection
TEST_ACTIVE = True
optical_encoder.start()
analog_in.start()
_thread.start_new_thread(overvoltage_protect_thread,())
if MODE_LIVE_READ:
    print_live_sensor_readout()
elif TIMED_DAQ_END == True:
    print(f'NIC Update:  the total duration of the test will be {time_DAQ_record} s.')

pulse_gen.start() # Start pulse generator last to trigger start of all other syncronized tasks
print("NIC Update: Data aquisition started.")
time.sleep(WAIT_TIME_PREINJECTION)
if TEST_ACTIVE: pump.infuse_const_rate(VOLUME_INJECTION, FLOW_RATE)
start_time_injection = time.perf_counter()
print("NIC Update: Pump successfully started.")

## Setting Method for Thread Termination
if TIMED_DAQ_END:
    while TEST_ACTIVE and time.perf_counter() - start_time_injection < (duration_injection + WAIT_TIME_POSTINJECTION):
        time.sleep(0.001)
    TEST_ACTIVE = False
else: 
    manual_stop_str = input("Press any key + Enter to stop the experiment: ")
    if manual_stop_str:
        TEST_ACTIVE = False
pulse_gen.stop()
print("NIC Update: Injection and data aquisition completed. Starting to clean-up DAQ tasks.")
optical_encoder.read_all(); analog_in.read_all() # Only necessary if using TDMS file reading with LoggingMode.LOG_AND_READ
optical_encoder.stop(); analog_in.stop()
pulse_gen.close(); optical_encoder.close(); analog_in.close()

print("NIC Update: Reading data from binary TDMS files.")
vol_data = TDMSHandler.read_tdms_data(optical_encoder) # (nL)
analog_in_data = TDMSHandler.read_tdms_data(analog_in) # (V)
pressure_data = analog_in_data[0] # (kPa)
force_data = analog_in_data[1] # (N)
if not len(vol_data) == len(pressure_data) == len(force_data):
    print("WARNING: Data lengths are not equal. Something has gone wrong.")
t = np.arange(len(vol_data))/SAMPLE_RATE

#### TEST OVERVIEW AND SAVIING ####
## Print Brief Test Overview
print("-----------------TEST OVERVIEW-----------------")
print(f'The maximum pressure was {max(pressure_data)} kPa')
print(f'The volume at maximum pressure was {vol_data[np.argmax(pressure_data)]} nL')
print("-----------------------------------------------")

## Write Data to a .csv File:
data_name = ["t","volume","pressure","force"]
data_unit = ["s","nL","kPa","N"]
data = [t, vol_data, pressure_data, force_data]
max_length = max(len(lst) for lst in data)
padded_data = [lst + [None] * (max_length - len(lst)) if max_length - len(lst) > 0 else lst for lst in data]
file = open(FILEPATH_CSV, 'w', newline ='')
with file:
    write = csv.writer(file)
    write.writerow(data_name)
    write.writerow(data_unit)
    write.writerows(zip(*padded_data))
file.close()
print(f"NIC Update: Data successfully written to {FILEPATH_CSV}.")