''' Volume-Controlled Cavitation Data Collection System (VCC-DCS)
Authors: Joseph Beckett, Brandon Brozich, and Esben Nielsen 
Description: Script for running needle-induced cavitation experiments in a custom experimental set-up
             developed at the University of Michigan under Prof. Jon Estrada starting in 2024. The script
             manages test definition, data acquisition, and data export to .csv files for additional 
             post-processing and inverse calibration of hyperelastic constituitive parameters in MATLAB.'''


#### INPUTS SECTION ####
## Test Method: Select the mode of operation for the script.
# LIVE_READ: Readout of sensor data is provided in real-time in the terminal. The pump is not actuated by the script. 
# VCC_CONSTANT_RATE: The pump is actuated to infuse a constant volume at a constant rate. Sensor data is recorded and saved to a .csv file.
# CALIBRATE_FRICTION: Run VCC_CONSTANT_RATE test with 5 repeated injections of 0.5 second durations with short pauses between to calibrate frictional forces for current FLOW_RATE.
MODE = "VCC_CONSTANT_RATE"
REPROCESS = False
APPLY_CORRECTIONS = False # If True, the csv will also include a corrected pressure reading without frictional additions.
RECORD_FORCE = False # If False, no force readings are taken

## System Defintions: Define parameters that only change with alterations to the physical equiptment.
DIAMETER_SYRINGE = 1.03     # Inner diameter of syringe (Hamilton Model 1705) in mm.
GAIN_PRESSURE_AMP = 262.78  # Should always be 368.04 unless DIP switches in the load cell Futek amplifier are changed from [0 0 0 0 0 0 1 0]
GAIN_FORCE_AMP = 262.78    # Should always be 1015.66 unless DIP switches in the load cell Futek amplifier are changed from [0 0 0 0 0 0 1 0]
VOLTAGE_RANGE_PRESSURE = [-10.0,10.0] # Voltage range (V) of pressure readings entering DAQ. Options limited to ±0.2 V, ±1 V, ±5 V, ±10 V.
VOLTAGE_RANGE_FORCE = [-10.0,10.0]   
LOAD_RANGE_PRESSURE = [-68.9,413.6] # (kPa) ~ (-10 to 60 psi) Pressure extrema for Pendotech PRESS-S-000.
LOAD_RANGE_FORCE = [-222.4,222.4] # (N) ~ (±50 lb) Force extrema for Futek LCM300.

## Experiment Settings
VOLUME_INJECTION = 100    # (nL) Amount of volume to be injected.
FLOW_RATE = 10             # (nL/s) Constant flow rate.

# Data Aquisition Settings
# Note: If analog input (i.e., force and pressure) readings are collected in excess of 1 kHz, ensure the 'Bandwidth' setting in the amplifier is set accordingly.
SAMPLE_RATE = 100            # (Hz) Sampling rate of all channels.                     
WAIT_TIME_PREINJECTION = 2  # (s) Duration that force and pressure readings are collected prior to starting the injection.
WAIT_TIME_POSTINJECTION = 2 # (s) Duration that force and pressure readings are collected after the injection is completed.
TIMED_DAQ_END = True    # If True, data aquisition ends after WAIT_TIME_POSTINJECTION, otherwise type a character in terminal and hit ENTER to manually 
                        # stop at any time after start of injection
from datetime import datetime
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M%S")  # Format: YYYYMMDD_HHMMSS
SPECIMENNAME = 'testing_new_motor'+TIMESTAMP
FOLDERPATH = 'C:\\Users\\jbecktt\\University of Michigan Dropbox\\Joseph Beckett\\um-research\\nic\\data\\troubleshooting\\'+SPECIMENNAME+'\\'
FOLDERPATH_CAL = 'C:\\Users\\jbecktt\\Documents\\nic-run-code\\cal\\'
FILEPATH_CSV = FOLDERPATH+SPECIMENNAME+'.csv'
CAL_FILE = FOLDERPATH_CAL + 'sfoil_friction_calibration_1_nLs-1_20241217_102956.csv'

#### CLASS DEFINITIONS ####
import serial, time, math, csv, nidaqmx, serial.tools.list_ports, os, sys, numpy as np, _thread, matplotlib.pyplot as plt, pandas as pd, keyboard
from nidaqmx.constants import AcquisitionType, TerminalConfiguration, EncoderType, AngleUnits, Edge, Level, FrequencyUnits, LoggingOperation, LoggingMode, Signal, READ_ALL_AVAILABLE
from nidaqmx.types import CtrFreq
from nptdms import TdmsFile
from itertools import chain
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
from collections import deque

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
                    if not i == meas_obj.data_type.index("force"):
                        voltage = tdms_file["analog_in"][meas_obj.ai_chan[i]].data
                        if raw == True:
                            data[i] = voltage #(V)
                        else:
                            if meas_obj.data_type[i] == "pressure":
                                data[i] = meas_obj._voltageV_to_pressurekPa(voltage) # (kPa)
                            elif meas_obj.data_type[i] == "force":
                                data[i] = meas_obj._voltageV_to_forceN(voltage) # (N)
                            else: raise Exception("Data type not supported yet for .tdms file reading for analog inputs.")
                return data
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
    def __init__(self, ci_chan, decoding_type, diameter_syringe, filepath, pulse_gen=None, sample_rate=None, n_samples=None, sample_mode=AcquisitionType.CONTINUOUS, savedata = True):
        super().__init__()
        units = AngleUnits.DEGREES
        pulses_per_rev = 5000 # E4T encoder spec
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
        if (isinstance(pulse_gen,PulseGenerator) or sample_rate is not None) and savedata == True:
            TDMSHandler.chkfile(filepath)
            self.task.in_stream.configure_logging(filepath,
                                                  LoggingMode.LOG_AND_READ,    
                                                  group_name='angular_disp',
                                                  operation=LoggingOperation.CREATE) # Create a new TDMS file. If the file already exists, NI-DAQmx returns an error.
        self.ci_chan = ci_chan
        self.diameter_syringe = diameter_syringe
        self.filepath = filepath
        if decoding_type == EncoderType.X_1:
            ppr_scale = 1.0
        elif decoding_type == EncoderType.X_2:
            ppr_scale = 2.0
        elif decoding_type == EncoderType.X_4:
            ppr_scale = 4.0
        else: raise Exception("Decoding type not recognized.")
        self.resolution = self._angle_deg_to_vol_nL(360.0/(ppr_scale*pulses_per_rev)) # nL/°

    def _angle_deg_to_vol_nL(self,angle_deg):
        #gear_ratio = 1/2.4
        #factor = gear_ratio*(1/24)*(25400/1)*(1/360)*10**-3 # mm/°of spindle rotation
        factor = 5*(1/360) # mm/°of spindle rotation
        A_syringe = math.pi*(self.diameter_syringe/2)**2 # mm^2
        vol = angle_deg*factor*A_syringe # mm^3 or μL
        return vol*10**3 # nL

    def read_volume(self,raw = False,all = False):
        cnt = self.read_all() if all else self.read()
        return self._angle_deg_to_vol_nL(np.array(cnt)).tolist() if not raw else cnt

class AnalogInput(NIDAQTask):
    def __init__(self, ai_chan, data_type, gain, ai_range_V, ai_range_conv, filepath, pulse_gen=None, sample_rate=None, n_samples=None, sample_mode=AcquisitionType.CONTINUOUS, savedata=True, force_bool = True):
        super().__init__() # Ensures that NIDAQTask Class is also initialized
        for i in range(len(ai_chan)):
            if force_bool or not data_type[i] == "force":
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
        if (isinstance(pulse_gen,PulseGenerator) or sample_rate is not None) and savedata == True:
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
    
    def read_pressure(self,raw=False,all=False):
        voltage = self.read_all() if all else self.read()
        return self._voltageV_to_pressurekPa(np.array(voltage[self.data_type.index("pressure")])).tolist() if not raw else voltage[self.data_type.index("pressure")]

    def _voltageV_to_forceN(self,voltage):
        N_lb = 4.44822 # N/lb conversion
        FSR = 50 # full-scale-range on loadcell (lb)
        Vexcit = 10 # excitation voltage of loadcell (V)
        RO = 1 # rated output of load cell (mV/V)
        sensitivity = FSR/(Vexcit*RO*10**-3*self.gain[self.data_type.index("force")]) # lb/V
        force = sensitivity*voltage*N_lb # N
        return force
    
    def read_force(self,raw=False,all=False):
        voltage = self.read_all() if all else self.read()
        return self._voltageV_to_forceN(np.array(voltage[self.data_type.index("force")])).tolist() if not raw else voltage[self.data_type.index("force")]

    def _get_overvoltage_limits(self):
        self.overvoltage_range = {}; self.overload = {}; self.overvoltage_range_conv = {}
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
            if self.data_type[i] == "force":
                self.overvoltage_range_conv[i] = self._voltageV_to_forceN(np.array(self.overvoltage_range[i])).tolist()
            elif self.data_type[i] == "pressure":
                self.overvoltage_range_conv[i] = self._voltageV_to_pressurekPa(np.array(self.overvoltage_range[i])).tolist()

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
        self.flow_rate_min = 0.069/27.5*(math.pi*(diameter/2)**2) # nL/s 
        self.flow_rate_max = 0.069/(26*10**-6)*(math.pi*(diameter/2)**2) # nL/s

        start_time_volume = time.perf_counter()
        self.query('svolume ' + '0' + ' ul') # Setting units to ul
        while not math.isclose(float(self.query('svolume').split()[0]), volume):
            self.query('svolume ' + str(volume) + ' ul')
            time.sleep(0.01)
            if time.perf_counter() - start_time_volume > 5: # 5 second timeout
                raise Exception('Unable to set syringe volume. Timed out')
            
    def __del__(self):
        if self.serial_port.is_open: self.serial_port.close()

class StepperMotor(NIDAQTask): 
    def __init__(self, port, ena_chan, diameter_syringe):
        super().__init__()
        try:
            self.serial_port = serial.Serial(port, 115200, timeout=1)
        except serial.SerialException as e:
            raise Exception(f"Failed to open serial port {port}: {e}")
        self.diameter_syringe = diameter_syringe
        self.task.do_channels.add_do_chan(ena_chan)
        self.task.write(False) # enable the motor

    def write(self, command):
        self.serial_port.write((command).encode())
    
    def infuse_const_rate(self):
        #self.write(f'RATE:{rate} VOLUME:{volume}')
        self.write('start')

    def disable_motor(self):
        self.task.write(True) # disable the motor
        self.write('stop')

#### TEST SET-UP ####
if MODE == "LIVE_READ":
    SAVE_DATA = False
    SAMPLE_RATE = 10e3
else: 
    TDMSHandler(FOLDERPATH)
    SAVE_DATA = True if REPROCESS == False else False
    if MODE == "CALIBRATE_FRICTION": SAMPLE_RATE = 5*FLOW_RATE

global analog_in, TEST_ACTIVE
## Initializing Syringe Pump and Sensors
if not MODE == "LIVE_READ":
    global stepper
    stepper = StepperMotor('COM7','/Dev2/port0/line0',DIAMETER_SYRINGE)
##pump = SyringePump('COM4')
##pump.set_syringe_dims(DIAMETER_SYRINGE, CAPACITY_SYRINGE)
duration_injection = VOLUME_INJECTION/FLOW_RATE
if MODE == "CALIBRATE_FRICTION" or MODE == "LIVE_READ":
    time_DAQ_record = 60 # Place holder (very long for now)
else:
    time_DAQ_record = WAIT_TIME_PREINJECTION + duration_injection + WAIT_TIME_POSTINJECTION
n_samples = math.ceil(time_DAQ_record*SAMPLE_RATE)
pulse_gen = PulseGenerator('/Dev2/PFI12','/Dev2/ctr1',SAMPLE_RATE,n_samples)
encoder_type = EncoderType.X_4 # X_2 is generally has lower resolution than X_4, but is less sensitive to vibrations
optical_encoder = OpticalEncoder('Dev2/ctr0', encoder_type, DIAMETER_SYRINGE, FOLDERPATH+'angular_disp.tdms', pulse_gen, savedata = SAVE_DATA)
analog_in = AnalogInput(['Dev2/ai0','Dev2/ai1'], ["pressure","force"], [GAIN_PRESSURE_AMP,GAIN_FORCE_AMP], 
                        [VOLTAGE_RANGE_PRESSURE,VOLTAGE_RANGE_FORCE], [LOAD_RANGE_PRESSURE,LOAD_RANGE_FORCE], 
                        FOLDERPATH+'analog_in.tdms', pulse_gen, savedata = SAVE_DATA, force_bool = RECORD_FORCE)

# Start system protection thread to protect from overvoltage of DAQ or overload of sensor(s)
def movmean(data, window_size):
    return pd.Series(data).rolling(window=window_size, min_periods=1).mean().to_numpy()
def overvoltage_protect_thread(live_read = False, refresh_rate = 10, record_force = RECORD_FORCE):
    global analog_in, TEST_ACTIVE
    if not live_read: global stepper
    if live_read:
        maxlen = int(np.ceil(SAMPLE_RATE/refresh_rate))
        data_queue_pressure = deque(maxlen=maxlen)
        if record_force: data_queue_force = deque(maxlen=maxlen)
    try:
        while TEST_ACTIVE:
            ai_data = analog_in.read_all()
            ctr_data = optical_encoder.read_volume(raw=False,all=True) # JGB: temporarily changed to raw=True
            if len(ai_data)> 0:
                safe_voltage = analog_in.safe_voltage_check(ai_data)
                if not safe_voltage:
                    stepper.stop_motor()
                    TEST_ACTIVE = False
                    raise Exception("Overvoltage protection triggered: analog voltage reading(s) exceeded safe range.")
            if MODE == "LIVE_READ" and TEST_ACTIVE:
                time.sleep(1/refresh_rate)
                vol_str = "Volume: %.3f nL  " % ctr_data[-1]
                if record_force:
                    data_queue_pressure.extend(ai_data[analog_in.data_type.index("pressure")])
                    data_queue_force.extend(ai_data[analog_in.data_type.index("force")])
                    force_V = np.mean(data_queue_force)
                    force_str = "Force: %.4f N (%.5f V)  " % (analog_in._voltageV_to_forceN(force_V), force_V)
                else: 
                    data_queue_pressure.extend(ai_data)
                pressure_V = np.mean(data_queue_pressure)
                pressure_str = "Pressure: %.4f kPa (%.5f V)  " % (analog_in._voltageV_to_pressurekPa(pressure_V), pressure_V)
                if record_force: 
                    print(f"{vol_str}\t{pressure_str}\t{force_str}", end="\r")
                else: 
                    print(f"{vol_str}\t{pressure_str}", end="\r")
    except KeyboardInterrupt:
        raise Exception("Live reading mode stopped by user.")
    #except Exception as e: 
    #    print(f"An error occurred during the overvoltage protection thread: {e}")
def write_csv(filename, data_name, data_unit, data):
    max_length = max(len(lst) for lst in data)
    padded_data = [lst + [None] * (max_length - len(lst)) if max_length - len(lst) > 0 else lst for lst in data]
    file = open(filename, 'w', newline ='')
    with file:
        write = csv.writer(file)
        write.writerow(data_name)
        write.writerow(data_unit)
        write.writerows(zip(*padded_data))
    file.close()
def read_csv(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        data_name = next(reader)
        data_unit = next(reader)
        data = [[] for _ in data_name]
        for row in reader:
            for i, value in enumerate(row):
                data[i].append(float(value) if value else None)
                
    return data_name, data_unit, data
def make_plot(filename, data_name, data_unit, data, volumeinjection = VOLUME_INJECTION, filt_window = None, title = None, record_force = RECORD_FORCE):
    rows = 3 if record_force else 2
    fig, axs = plt.subplots(rows, 1, figsize=(10, 12))
    axs[0].plot(data[data_name.index("time")], data[data_name.index("volume")])
    axs[0].set_ylabel('Volume ('+data_unit[data_name.index("volume")]+')')
    axs[0].set_xlabel('Time ('+data_unit[data_name.index("time")]+')')
    axs[1].plot(data[data_name.index("time")], data[data_name.index("pressure")],label='Raw Data')
    if filt_window is not None:
        for i in range(len(filt_window)):
            if round(SAMPLE_RATE*volumeinjection/FLOW_RATE*filt_window[i]) > 2:
                axs[1].plot(data[data_name.index("time")], savgol_filter(data[data_name.index("pressure")],round(SAMPLE_RATE*volumeinjection/FLOW_RATE*filt_window[i]),2),label=f'SG Window: {filt_window[i]} s')
                if APPLY_CORRECTIONS == True and "pressure_corrected" in data_name:
                    axs[1].plot(data[data_name.index("time")], savgol_filter(data[data_name.index("pressure_corrected")],round(SAMPLE_RATE*volumeinjection/FLOW_RATE*filt_window[i]),2),label=f'Corrected Data. SG: {filt_window[i]} s')
    axs[1].set_ylabel('Pressure ('+data_unit[data_name.index("pressure")]+')')
    axs[1].legend()
    axs[1].set_xlabel('Time ('+data_unit[data_name.index("time")]+')')
    if record_force:
        axs[2].plot(data[data_name.index("time")], data[data_name.index("force")])
        axs[2].set_ylabel('Force '+data_unit[data_name.index("force")])
        axs[2].set_xlabel('Time ('+data_unit[data_name.index("time")]+')')
    if title is not None: plt.suptitle(title, fontsize=16, fontweight='bold')
    plt.savefig(filename)
    plt.show()
    if APPLY_CORRECTIONS == True and MODE == "CALIBRATE_FRICTION":
        residual = np.mean(data[data_name.index("pressure_corrected")])
        print(f"Mean of pressure_corrected is: {residual}")
def friction_calibration(data_name, data, N_CYCLES, T_INJECT, T_PAUSE_CYCLE):
    cnt_pace = math.ceil(optical_encoder.resolution*SAMPLE_RATE/FLOW_RATE) # S/cnt
    mm_fraction = 0.05
    window_mm = round(mm_fraction*T_INJECT*SAMPLE_RATE) # Adjust window size based on total samples collected in injection time -- not completely right
    # N_0 = round(0.25*T_INJECT*SAMPLE_RATE)
    N_0_wait = round(0.05*T_PAUSE_CYCLE*SAMPLE_RATE)
    N_0_bigwait = round(0.50*T_PAUSE_CYCLE*SAMPLE_RATE)
    N_0_inject = round(0.25*T_INJECT*SAMPLE_RATE)

    idx_time = data_name.index("time"); idx_vol = data_name.index("volume"); idx_pressure = data_name.index("pressure"); idx_force = data_name.index("force")
    vdot_raw = np.diff(data[idx_vol])/np.diff(data[idx_time])
    vdot = savgol_filter(vdot_raw,window_mm,2)
    v2dot = savgol_filter(np.diff(vdot)/np.diff(data[idx_time][1:]),window_mm,2)
    vdot_is0 = vdot == 0
    vdot_raw_is0 = vdot_raw == 0
    vdot_raw_lt0 = vdot_raw < 0
    vdot_raw_gt0 = vdot_raw > 0
    vdot_raw_le0 = vdot_raw <= 0

    # Find start of each injection cycle
    LIMITS_FOUND = False; itr = 0; N = N_0_wait
    while LIMITS_FOUND == False:
        idx_starts = [] # index of the start of each injection cycle
        itr += 1
        for i in range(N,len(vdot)-N):
            if vdot_raw_gt0[i] == True and all(vdot_raw_is0[i-N:i]):
                idx_starts.append(i)
        if len(idx_starts) != N_CYCLES:
            if itr == 30 or N == 1:
                raise Exception("Friction calibration failed: Could not find bounds of the 5 injections.")
            elif len(idx_starts) > N_CYCLES: 
                N = round(1.1*N)
            else: N = round(N*0.5)
        else: 
            LIMITS_FOUND = True
            print(f"Found idx_starts in {itr} iterations. N/N_0_wait = {N/N_0_wait}")
    def isNextNonzeroValueNegative(vec):
        idx_first_nonzero = np.where(np.array(vec) != 0)[0][0]
        if vec[idx_first_nonzero] < 0:
            return True
        else: return False

    # Find end of each injection cycle
    LIMITS_FOUND = False; itr = 0; N = N_0_wait
    while LIMITS_FOUND == False:
        idx_stops = [] # index of the start of each injection cycle
        itr += 1
        for i in range(N,len(vdot)-N+1):
            if vdot_raw_gt0[i-1] == True and any(vdot_raw_lt0[i-N:i-1]) == False and (sum(vdot_raw_gt0[i-10*cnt_pace:i])>(10-1)) and (all(vdot_raw_le0[i:i+N]) or isNextNonzeroValueNegative(vdot_raw[i:i+N])): #and vdot_raw_is0[i-1-cnt_pace] == False:
            # if vdot_raw_is0[i-1] == False and np.mean(vdot_raw_le0[i-cnt_pace:i]):
                idx_stops.append(i)
        if len(idx_stops) != N_CYCLES:
            if itr == 30 or N == 1:
                raise Exception("Friction calibration failed: Could not find bounds of the 5 injections.")
            elif len(idx_stops) > N_CYCLES: 
                N = math.ceil(1.1*N)
            else: N = math.floor(0.5*N)
        else: 
            LIMITS_FOUND = True
            print(f"Found idx_stops in {itr} iterations. N/N_0_wait = {N/N_0_wait}")
    # Find start of each steady state region
    LIMITS_FOUND = False; itr = 0; C = 1
    while LIMITS_FOUND == False:
        idx_ss_starts = [] # index of the start of each steady state region
        itr += 1
        speed_up = v2dot > np.std(v2dot)*C
        for i in range(N_0_bigwait,len(vdot)-N_0_inject):
            if speed_up[i-1] == True and all(speed_up[i:i+N_0_inject] == False) and all(vdot_is0[i:i+N_0_inject] == False):
                idx_ss_starts.append(i)
        if len(idx_ss_starts) != N_CYCLES:
            if itr == 30:
                raise Exception("Friction calibration failed: Could not find bounds of the 5 steady state regions.")
            elif len(idx_ss_starts) > N_CYCLES:
                C *= 1.1
            elif len(idx_ss_starts) < N_CYCLES:
                C *= 0.9
        else: 
            LIMITS_FOUND = True
            print(f"Found idx_ss_starts in {itr} iterations.")
    # Find end of each steady state region
    LIMITS_FOUND = False; itr = 0; C = 1
    while LIMITS_FOUND == False:
        idx_ss_stops = [] # index of the end of each steady state region
        itr += 1
        slow_down = v2dot < -np.std(v2dot)*C
        for i in range(N_0_inject,len(vdot)-N_0_bigwait):
            if slow_down[i+1] == True and all(slow_down[i-N_0_inject+1:i+1] == False):
                idx_ss_stops.append(i)
        if len(idx_ss_stops) != N_CYCLES:
            if itr == 20:
                raise Exception("Friction calibration failed: Could not find bounds of the 5 steady state regions.")
            else: C *= 1.1
        else: 
            LIMITS_FOUND = True
            print(f"Found idx_ss_stops in {itr} iterations.")

    # Segement data into 5 individual injections
    data_segment = []
    for start, stop in zip(idx_starts, idx_stops):
        data_segment.append([data[i][start:stop] for i in range(len(data))])
    pressure_baseline = np.mean(data[idx_pressure][0:idx_starts[0]])
    pressure_baseline_std = np.std(data[idx_pressure][0:idx_starts[0]])
    # force_baseline = np.mean(data[idx_force][0:idx_starts[0]])
    # force_baseline_std = np.std(data[idx_force][0:idx_starts[0]])

    idx_ss_starts_segment = []; idx_ss_stops_segment = []
    for i in range(N_CYCLES):
        data_segment[i][idx_time] -= data[idx_time][idx_starts[i]-1] # correct for time offset between segments
        data_segment[i][idx_vol] -= data[idx_vol][idx_starts[i]-1] # correct for volume offset between segments
        data_segment[i][idx_pressure] -= pressure_baseline # correct for pressure baseline
        idx_ss_starts_segment.append(idx_ss_starts[i]-idx_starts[i]) # cropped to a conservative start-up region
        idx_ss_stops_segment.append(idx_ss_stops[i]-idx_starts[i]) # cropped to a conservative start-up region

    segment_lens = [len(data_segment[i][0]) for i in range(N_CYCLES)]
    data_avg = [[] for _ in range(len(data))]
    pressure_filt = [savgol_filter(data_segment[i][idx_pressure],window_mm,2) for i in range(N_CYCLES)]
    # force_filt = [savgol_filter(data_segment[i][idx_force],window_mm,2) for i in range(N_CYCLES)]
    for i in range(min(segment_lens)):
        #fit_avg.append(np.mean([fit[k][i] for k in range(N_CYCLES)]))
        for j in range(len(data)):
            data_avg[j].append(np.mean([data_segment[k][j][i] for k in range(N_CYCLES)]))
    pressure_avg_filt = savgol_filter(data_avg[idx_pressure],window_mm,2)
    # force_avg_filt = savgol_filter(data_avg[idx_force],window_mm,2)

    ss_start_est = max(idx_ss_starts_segment)
    ss_stop_est = min(idx_ss_stops_segment)
    ss_safe_len = ss_stop_est - ss_start_est + 1
    ss_start_safe = ss_start_est + round(0.05*ss_safe_len)
    ss_stop_safe = ss_stop_est - round(0.05*ss_safe_len)

    pressure_drop_ss = np.mean(pressure_avg_filt[ss_start_safe:ss_stop_safe])

    residual_ss = [np.array(pressure_filt[i][ss_start_safe:ss_stop_safe] - pressure_drop_ss) for i in range(len(pressure_filt))]
    residual_trans = [np.array(pressure_filt[i][0:ss_start_safe]-pressure_avg_filt[0:ss_start_safe]) for i in range(len(pressure_filt))]
    mean_residual_ss = float(np.mean(np.concatenate(residual_ss)))
    mean_residual_trans = float(np.mean(np.concatenate(residual_trans)))
    std_residual_ss = float(np.std(np.concatenate(residual_ss)))
    std_residual_trans = float(np.std(np.concatenate(residual_trans)))

    print(f"Mean ± SD pressure residual (ss//trans): {mean_residual_ss} ± {std_residual_ss} kPa // {mean_residual_trans} ± {std_residual_trans} kPa")
    
    # Make a plot of the data_segments but they're all on the same plot (not a subplot)
    # Make all lines shades of blue
    #plt.figure(figsize=(10,6))
    fig, axs = plt.subplots(3, 2, figsize=(20, 20))
    plt.suptitle(f'Frictional Pressure Drop Calibration. Flow Rate: {FLOW_RATE} nL/s', fontsize=16, fontweight='bold')
    for i in range(len(data_segment)):
        axs[0,0].plot(data_segment[i][idx_time],data_segment[i][idx_vol],label=f'Injection {i+1}',color=(0,0,1,(i+1)/(len(data_segment)+1)))
        axs[1,0].plot(data_segment[i][idx_time],vdot[idx_starts[i]:idx_stops[i]],label=f'Injection {i+1}',color=(0,0,1,(i+1)/(len(data_segment)+1)))
        axs[2,0].plot(data_segment[i][idx_time],data_segment[i][idx_pressure],label=f'Injection {i+1}',color=(0,0,1,(i+1)/(len(data_segment)+1)))
        axs[0,1].plot(data_segment[i][idx_time],pressure_filt[i],label=f'Injection {i+1}',color=(0,0,1,(i+1)/(len(data_segment)+1)))
        axs[2,1].plot(data_avg[idx_time][ss_start_safe:ss_stop_safe],residual_ss[i],color=(0,0,1,(i+1)/(len(data_segment)+1)))
        axs[2,1].plot(data_avg[idx_time][0:ss_start_safe],residual_trans[i],color=(0.3,0.6,1,(i+1)/(len(data_segment)+1)))
    axs[1,1].plot(data_avg[idx_time],pressure_avg_filt,label=f'Injection {i+1}',color=(0,0,1,(i+1)/(len(data_segment)+1)))
    axs[1,1].axvline(data_avg[idx_time][max(idx_ss_starts_segment)],color=(0,0,1,0.5),linestyle='--')
    axs[1,1].axvline(data_avg[idx_time][min(idx_ss_stops_segment)],linestyle='--',color=(0,0,1,0.5))
    axs[1,1].axvline(data_avg[idx_time][ss_start_safe],linestyle='-',color=(0,0,1,0.5))
    axs[1,1].axvline(data_avg[idx_time][ss_stop_safe],linestyle='-',color=(0,0,1,0.5))
    axs[0,0].set_ylabel('Volume (nL)'); axs[0,0].set_xlabel('Time (s)')
    axs[1,0].set_ylabel('Flow Rate (nL/s)'); axs[1,0].set_xlabel('Time (s)')
    axs[2,0].set_ylabel('Pressure (kPa)'); axs[2,0].set_xlabel('Time (s)'); axs[2,0].legend()
    axs[0,1].set_ylabel('Pressure (kPa) (SG filter)'); axs[0,1].set_xlabel('Time (s)')
    axs[1,1].set_ylabel('Avg Pressure (kPa) (SG filter)'); axs[1,1].set_xlabel('Time (s)'); axs[1,1].legend()
    axs[2,1].set_ylabel('Pressure Residual (kPa)'); axs[2,1].set_xlabel('Time (s)')
   
    plt.savefig(FOLDERPATH+f'CVE_Frictional_Calibration_{FLOW_RATE}_nLs-1.png')
    plt.show()

    if abs(np.mean(pressure_avg_filt)) < 0.001: # and np.std(savgol_filter(data_avg[idx_pressure],window_mm,2)) < 0.002:
        PRESSURE_DROP_CORRECTION = [0] # Pressure drop correction is not necessary because noise exceeds pressure drop
        t_ss_start = [0]
        pressure_drop_ss = [0]
        t_trans = []
        pressure_trans = []
        # force_trans = []
    else: 
        print("Automated calibration has completed. It believes there is meaningful pressure drop. Press 1 and Enter to save the calibration.")
        PRESSURE_DROP_CORRECTION = [1]
        t_ss_start = [data_avg[idx_time][ss_start_safe]] # beyond this time the steady state pressure drop can be used to correct the data
        pressure_drop_ss = [pressure_drop_ss] # already defined above plotting
        t_trans = data_avg[idx_time][0:ss_start_safe]
        pressure_trans = pressure_avg_filt[0:ss_start_safe]
        # force_trans = force_avg_filt[0:ss_start_safe]
    if PRESSURE_DROP_CORRECTION == [0]:
        msg = "The code believes that the noise in the pressure sensor exceeds the pressure drop. No correction will be applied."
    else: msg = f"The code believes that a pressure drop correction is necessary at {FLOW_RATE} nL/s."
    user_input = input("Automated calibration has completed. "+msg+" Press 1 and hit Enter to save the calibration: ")
    if user_input == "1":
        data_name_cal = ["PRESSURE_DROP_CORRECTION","t_ss_start","pressure_drop_ss","t_trans","pressure_trans","mean_residual_ss","mean_residual_trans","std_residual_ss","std_residual_trans"]
        data_unit_cal = ["bool","s","kPa","s","kPa","kPa","kPa","kPa","kPa"]
        data_cal = [PRESSURE_DROP_CORRECTION, t_ss_start, pressure_drop_ss, t_trans, pressure_trans, [mean_residual_ss], [mean_residual_trans], [std_residual_ss], [std_residual_trans]]
        write_csv(FOLDERPATH+f'friction_calibration_{FLOW_RATE}.csv', data_name_cal, data_unit_cal, data_cal)
        write_csv(FOLDERPATH_CAL+f'friction_calibration_{FLOW_RATE}_nLs-1_{TIMESTAMP}.csv', data_name_cal, data_unit_cal, data_cal)
        print("Friction calibration data written to csv file.")
    else: print("Calibration aborted. No calibration file saved.")
    return PRESSURE_DROP_CORRECTION
def friction_correction(t,volume,pressure,cal_file):
    pressure_correct = pressure.copy()
    cal_name, cal_unit, cal_data = read_csv(cal_file)
    t_ss_start = cal_data[cal_name.index("t_ss_start")][0] # This is time *after* start of injection
    idx_inject_start = find_injection_start(volume)
    idx_inject_end = find_injection_end(volume)
    t_inject_start = t[idx_inject_start]
    t_inject_end = t[idx_inject_end]
    idx_ss_start = np.where(t >= t_ss_start + t_inject_start)[0][0]
    idx_ss_end = np.where(t >= t_inject_end)[0][0]
    trans_interp = interp1d(cal_data[cal_name.index("t_trans")],cal_data[cal_name.index("pressure_trans")],fill_value="extrapolate")
    pressure_init = np.mean(pressure[:idx_inject_start])

    pressure_correct -= pressure_init
    pressure_correct[idx_inject_start:idx_ss_start] -= trans_interp(t[idx_inject_start:idx_ss_start]-t[idx_inject_start])
    pressure_correct[idx_ss_start:idx_ss_end] -= cal_data[cal_name.index("pressure_drop_ss")][0]
    return pressure_correct
def find_injection_start(volume, N_CYCLES=1):
    N_0_wait = round(0.05*WAIT_TIME_PREINJECTION*SAMPLE_RATE)
    vdot = np.diff(volume)
    vdot_is0 = vdot == 0
    vdot_gt0 = vdot > 0
    LIMITS_FOUND = False; itr = 0; N = N_0_wait
    while LIMITS_FOUND == False:
        idx_starts = [] # index of the start of each injection cycle
        itr += 1
        for i in range(N,len(vdot)-N):
            if vdot_gt0[i] == True and all(vdot_is0[i-N:i]):
                idx_starts.append(i)
        if len(idx_starts) != N_CYCLES:
            if itr == 30 or N == 1:
                raise Exception("Friction calibration failed: Could not find bounds of the injection(s).")
            elif len(idx_starts) > N_CYCLES: 
                N = round(1.5*N)
            else: N = round(N*0.5)
        else: 
            LIMITS_FOUND = True
            print(f"Found idx_starts in {itr} iterations. N/N_0_wait = {N/N_0_wait}")
            return idx_starts[0] if N_CYCLES == 1 else idx_starts
def find_injection_end(volume, N_CYCLES=1):
    cnt_pace = math.ceil(optical_encoder.resolution*SAMPLE_RATE/FLOW_RATE) # S/cnt
    def isNextNonzeroValueNegative(vec):
        idx_first_nonzero = np.where(np.array(vec) != 0)[0][0]
        if vec[idx_first_nonzero] < 0:
            return True
        else: return False
    vdot = np.diff(volume)
    vdot_lt0 = vdot < 0
    vdot_le0 = vdot <= 0
    vdot_gt0 = vdot > 0
    N_0_wait = round(0.05*WAIT_TIME_PREINJECTION*SAMPLE_RATE)
    # Find end of each injection cycle
    LIMITS_FOUND = False; itr = 0; N = N_0_wait
    while LIMITS_FOUND == False:
        idx_stops = [] # index of the start of each injection cycle
        itr += 1
        for i in range(N,len(vdot)-N+1):
            if vdot_gt0[i-1] == True and any(vdot_lt0[i-N:i-1]) == False and (sum(vdot_gt0[i-10*cnt_pace:i])>(10-1)) and (all(vdot_le0[i:i+N]) or isNextNonzeroValueNegative(vdot[i:i+N])): #and vdot_raw_is0[i-1-cnt_pace] == False:
            # if vdot_raw_is0[i-1] == False and np.mean(vdot_raw_le0[i-cnt_pace:i]):
                idx_stops.append(i)
        if len(idx_stops) != N_CYCLES:
            if itr == 30 or N <= 1:
                raise Exception("Friction calibration failed: Could not find bounds of the injection(s).")
            elif len(idx_stops) > N_CYCLES: 
                N = math.ceil(1.1*N)
            else: N = math.ceil(0.5*N)
        else: 
            LIMITS_FOUND = True
            print(f"Found idx_stops in {itr} iterations. N/N_0_wait = {N/N_0_wait}")
            return idx_stops[0] if N_CYCLES == 1 else idx_stops
            
class NICTestMethod:
    def __init__(self,test,reprocess=False):
        self.test = test
        if test == "LIVE_READ":
            self.print_live_sensor_readout()
        elif test == "VCC_CONSTANT_RATE":
            self.vcc_constant_rate(reprocess=reprocess)
        elif test == "CALIBRATE_FRICTION":
            self.vcc_constant_rate(calibrate_friction=True,reprocess=reprocess)
        
    def print_live_sensor_readout(self,refresh_rate=10):
        global TEST_ACTIVE
        optical_encoder.start(); analog_in.start(); pulse_gen.start()
        try: 
            _thread.start_new_thread(overvoltage_protect_thread(live_read=True,refresh_rate=refresh_rate),())
        except Exception as e:
            if TEST_ACTIVE: TEST_ACTIVE = False
            print(f"\n{e}")
        finally: 
            pulse_gen.close(); optical_encoder.close(); analog_in.close()

    def vcc_constant_rate(self,calibrate_friction=False,reprocess=False,record_force=RECORD_FORCE):
        if reprocess == False:
            global TEST_ACTIVE
            ## Prepare DAQ for test
            if TIMED_DAQ_END:
                print(f'NIC Update:  the total duration of the test will be {time_DAQ_record} s.')
            optical_encoder.start()
            analog_in.start()
            _thread.start_new_thread(overvoltage_protect_thread,())
            pulse_gen.start()
            ## Start of data aquisition and test
            #pulse_gen.start() # Start pulse generator last to trigger start of all other syncronized tasks
            print("NIC Update: Data aquisition started.")
            if not calibrate_friction:
                time.sleep(WAIT_TIME_PREINJECTION)
                #if TEST_ACTIVE: stepper.infuse_const_rate(VOLUME_INJECTION, FLOW_RATE)
                if TEST_ACTIVE: stepper.infuse_const_rate()
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
            elif calibrate_friction:
                N_CYCLE = 3
                T_INJECT = optical_encoder.resolution/FLOW_RATE*1e3*0.25 # Time for 1000 optical encoder pulses to occur
                VOLUME_INJECTION_CAL = FLOW_RATE*T_INJECT
                T_PAUSE_INIT = max([1 + 0.1*T_INJECT,2]) # Modified by T_INJECT to ensure no issues occur even if sampling rate is low for low rate tests
                T_PAUSE_CYCLE = T_PAUSE_INIT
                print(f"NIC Update: Calibration will take approx. {T_PAUSE_INIT+N_CYCLE*(T_PAUSE_CYCLE+T_INJECT)} s to complete.")
                t_startcal = time.perf_counter()
                n = 0 
                time.sleep(T_PAUSE_INIT)
                while n < N_CYCLE: # run 5 cycles
                    n += 1
                    t_startpump = time.perf_counter()
                    if TEST_ACTIVE: stepper.infuse_const_rate()
                    t_endpump = time.perf_counter()
                    time.sleep(T_INJECT-(t_endpump-t_startpump)+T_PAUSE_CYCLE) # effectively injecting and waiting T_PAUSE_CYCLE seconds before ending cycle
                TEST_ACTIVE = False
                t_endcal = time.perf_counter()
                print(f"NIC Update: Total time of test: {t_endcal-t_startcal} s.")
            pulse_gen.stop()
            print("NIC Update: Injection and data aquisition completed. Starting to clean-up DAQ tasks.")
            optical_encoder.read_all(); analog_in.read_all() # Only necessary if using TDMS file reading with LoggingMode.LOG_AND_READ
            optical_encoder.stop(); analog_in.stop()
            pulse_gen.close(); optical_encoder.close(); analog_in.close(); stepper.close()

        print("NIC Update: Reading data from binary TDMS files.")
        vol_data = TDMSHandler.read_tdms_data(optical_encoder) # (nL)
        analog_in_data = TDMSHandler.read_tdms_data(analog_in) # (V)
        pressure_data = analog_in_data[analog_in.data_type.index("pressure")] # (kPa)
        if record_force:
            force_data = analog_in_data[analog_in.data_type.index("force")] # (N)
            len_check = len(vol_data) == len(pressure_data) == len(force_data)
        else: 
            len_check = len(vol_data) == len(pressure_data)
        if not len_check:
            print("WARNING: Data lengths are not equal. Something has gone wrong.")
        t = np.arange(len(vol_data))/SAMPLE_RATE
        ## Print Brief Test Overview
        print("-----------------TEST OVERVIEW-----------------")
        print(f'The maximum pressure was {max(pressure_data)} kPa')
        print(f'The volume at maximum pressure was {vol_data[np.argmax(pressure_data)]} nL')
        print("-----------------------------------------------")
        ## Write Data to a .csv File:
        data_name = ["time", "volume", "pressure"]
        data_unit = ["s", "nL", "kPa"]
        data = [t, vol_data, pressure_data]
        if record_force:
            data_name.append("force")
            data_unit.append("N")
            data.append(force_data)
        if APPLY_CORRECTIONS and MODE != "CALIBRATE_FRICTION":
            try:
                pressure_data_corrected = friction_correction(t,vol_data,pressure_data,CAL_FILE)
                data_name.append("pressure_corrected")
                data_unit.append("kPa")
                data.append(pressure_data_corrected)
            except Exception as e: 
                print(f"An error occurred during friction calibratino: {e}")
        if calibrate_friction: 
            friction_calibration(data_name,data,N_CYCLE,T_INJECT,T_PAUSE_CYCLE)

        write_csv(FILEPATH_CSV, data_name, data_unit, data)
        print(f"NIC Update: Data successfully written to {FILEPATH_CSV}.")
        FILEPATH_PNG = FILEPATH_CSV.replace(".csv", ".png")
        if not calibrate_friction: 
            make_plot(FILEPATH_PNG, data_name, data_unit, data, filt_window = [0.025, 0.05], title= f"FLOW_RATE: {FLOW_RATE} nL/s. SAMPLE_RATE: {SAMPLE_RATE} Hz.")
            print(f"NIC Update: PLOT successfully written to {FILEPATH_PNG}.")

#### RUN VCC TEST ####
## Start Data Acquisition and Injection
TEST_ACTIVE = True
if MODE in ["LIVE_READ", "VCC_CONSTANT_RATE", "CALIBRATE_FRICTION"]:
    NICTestMethod(MODE,REPROCESS)
else: 
    raise Exception("Invalid test mode selected. Please select a valid test mode.")