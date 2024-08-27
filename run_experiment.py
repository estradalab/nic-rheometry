# Name: Volume-Controlled Cavitation Data Collection System (VCC-DCS)
# Author: Brandon Brozich
# Date: 7/17/2024
# Description:
#   Controls the volume-controlled cavitation devices for running experiments. Outputs volume readings (nL) in the first row of output.csv,
#   and outputs pressure readings (kPa) in the third row of output.csv.

TARE_SENSOR = False        # Set to "True" if you want the pressure readings to be zeroed upon pressing the Run button. Set to "False" otherwise.
LIVE_PRESSURE_READING_MODE = False    # Set to "True" if you would like to view live pressure readings without moving the syringe pump.
LIVE_FORCE_READING_MODE = True
                            # To stop receiving live-readings, click in the terminal where the readings are being output, and hit Ctrl-C.
                            # Set to "False" if you do not want to view live pressure readings and instead want to run the syringe pump.
SYRINGE_DIAMETER = 1.03     # Innder diameter of syringe in millimeters (mm).
SYRINGE_CAPACITY = 50       # Syringe volumetric capacity in microliters (uL).
TARGET_VOLUME = 5000       # Amount of volume to be injected in nanoliters (nL).
FLOW_RATE = 40              # Flow rate (for constant flow rate test) in nanoliters per second (nL/sec).
MAXIMUM_PRESSURE = 300      # Maximum pressure (kPa) before test is stopped, to prevent pressure sensor overpressuring.
                            # Warning: Pressures may exceed this number if the ADC's maximum operating range is exceeded.
ADC_GAIN = 32               # Gain of the ADC as set in the code.py file. Must be either 1, 2, 4, 8, 16, 32, 64, or 128.
FUTEK_AMP_GAIN = 1015.66 # Should always be 1015.66 unless DIP switches in amplifier are changed from [0 0 0 0 1 1 1 0]

#FILENAME = '20240808_NIC_Gelatin10percent-cuvette2_test_0flow.csv'
#FILENAME = '20240809_NIC_Gelatin10percent-cuvette7_test1_withRainX_1hr_15min_wait_2000nLps.csv'
#FILENAME = '20240809_compliancetestbubble_iinair2blahblah.csv'
#FILENAME = '20240809_NIC_Gelatin10percent_cuvette8_test1_curedaroundneedle_RainX_curedmaybe1_5hr_godeeperandgoagain'
#FILENAME = '20240812_NIC_Gelatin10percent_cuvetterando_test2_satinfixtureallweekend_probablydry.csv'
#FILENAME = '20240813_NIC_PA5percent0812_cuvette1_test1_APEScoating_18ishhrwait_neverfractured.csv'
#FILENAME = '20240819_NIC_PA5percent0816_cuvette14_test1_uncoated25Gsharp_60minwait_rate400_vol10k_retraction1mm.csv'
#FILENAME = '20240821_LoadCellTesting_test1.csv'
#FILENAME = 'watertestdumbgoodforce_10.csv'
FILENAME = '20240824_NIC_PA5percent0816_cuvette4_test1_uncoated25Gsharp_5minwait_rate40_vol5k_retraction2mm.csv'
#NOFLOWREAD = True

import serial
import time
import math
import csv
import sys
import nidaqmx # added for loadcell stuff
import concurrent.futures

class PressureSensor:
    def set_gain(self, gain):
        '''Sets gain of ADC (must be either 1, 2, 4, 8, 16, 32, 64, or 128).
            WARNING: This function is slow (~0.5 seconds to respond.)'''
        self.serial_port.write(('GAIN:' + str(gain) + '\r').encode())
        start = time.perf_counter()
        while time.perf_counter() - start < 5: # Time out after 5 seconds
            line = self.serial_port.readline().strip()
            if line == b'': # Skip if line is blank
                continue
            if line[-1] == 0x03:
                if line.strip(b'\x03')  == str(gain).encode():
                    time.sleep(0.1) # This is here because the ADC does not respond to gain changes immediately.
                    return
        raise Exception("Unable to set ADC gain.")

    def __init__(self, port):
        self.serial_port = serial.Serial(port, 115200, timeout=1)
        self.set_gain(ADC_GAIN)

    def convert_to_kpa(data_y, gain):
        mV_psi_V = 0.2584 # 0.2584 mV/psi/V
        V_ref = 3.0 # reference voltage
        fullScale_factor = 2
        mV_psi = mV_psi_V * V_ref # mV/psi
        V_psi = mV_psi / 1000 # V/psi
        V_kPa = V_psi / 6.89476 # V/kPa
        nBits = 2**24 # 24-bit ADC, i.e. 2^24 bits
        measuredVoltage = (data_y / nBits) * V_ref
        actualVoltage = (measuredVoltage  / gain) / fullScale_factor
        data_y_converted = (actualVoltage / V_kPa)
        return data_y_converted
    
    def get_pressure(self):
        '''Returns float of pressure reading in kilopascal (kPa).'''
        self.serial_port.write(b'READ\r')
        start = time.perf_counter()
        while time.perf_counter() - start < 5: # Time out after 5 seconds
            line = self.serial_port.readline().strip()
            if line == b'': # Skip if line is blank
                continue
            if line[-1] == 0x03:
                return PressureSensor.convert_to_kpa(int(line.strip(b'\x03')), ADC_GAIN)
        raise Exception("Unable to get pressure from sensor. Timed out.")

    def tare(self):
        '''Tares the pressure sensor.'''
        self.serial_port.write(b'TARE\r')
        start = time.perf_counter()
        while time.perf_counter() - start < 5: # Time out after 5 seconds
            line = self.serial_port.readline().strip()
            if line == b'': # Skip if line is blank
                continue
            if line[-1] == 0x03:
                if line.strip(b'\x03')  == b'SUCCESS':
                    return
        raise Exception("Unable to tare pressure sensor.")

    def __del__(self):
        if self.serial_port.is_open: self.serial_port.close()

# Load Cell Class (to be implemented if necessary):
    # function for getting force
    # function to tare
    # variable for storing the COM port
class NIDAQTask:
    def __init__(self):
        self.task = nidaqmx.Task()

    def __del__(self):
        self.task.close()

    def add_ai_channel(self, channel_handle, RSE_ON):
        if RSE_ON:
            tconfig = nidaqmx.constants.TerminalConfiguration.RSE
        else:
            tconfig = nidaqmx.constants.TerminalConfiguration.DEFAULT

        self.task.ai_channels.add_ai_voltage_chan(physical_channel=channel_handle, terminal_config=tconfig)

    def read(self):
        return self.task.read()
    
    def print_read(self):
        print(self.read())    

class LoadCell(NIDAQTask):
    def __init__(self, channel_handle, RSE_ON, gain):
        NIDAQTask.__init__(self)
        NIDAQTask.add_ai_channel(self, channel_handle, RSE_ON)
        self.gain = gain

    def convert_to_N(self,voltage):
        N_lb = 4.44822 # N/lb conversion
        FSR = 50 # full-scale-range on loadcell (lb)
        Vexcit = 10 # excitation voltage of loadcell (V)
        RO = 1 # rated output of load cell (mV/V)
        sensitivity = FSR/(Vexcit*RO*10**-3*self.gain) # lb/V
        force = sensitivity*voltage*N_lb # N
        return force

    def get_voltage(self):
        return self.read()
    
    def get_force(self):
        return self.convert_to_N(self.read())

    def print_voltage(self):
        print("%.3f V" % self.read())

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
        self.serial_port = serial.Serial(port, 115200, timeout=1)
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

# Initializing pump & pressure sensor:
pump = SyringePump('COM3') # Syringe pump has ~3.5 ms response time
pressure_sensor = PressureSensor('COM7') # Pressure sensor has ~12 ms response time (can be a little faster but becomes a bit unstable)
load_cell = LoadCell('Dev9/ai0', True, FUTEK_AMP_GAIN)

# Tare pressure readings:
if TARE_SENSOR:
    pressure_sensor.tare()
    print("Pressure Sensor tared.")
    sys.exit()

# While loop for seeing pressure readings live (hit Ctrl-C in the terminal to stop):
while LIVE_PRESSURE_READING_MODE:
    print(f"Pressure:\t {pressure_sensor.get_pressure()} kPa")
    time.sleep(0.25)

while LIVE_FORCE_READING_MODE:
    print(f"Force:\t {load_cell.get_force()} N")
    time.sleep(0.25)

volume_time = []
volume_data = []
pressure_time = []
pressure_data = []
force_time = []
force_data = []

pump.set_syringe_dims(SYRINGE_DIAMETER, SYRINGE_CAPACITY)
pump.infuse_const_rate(TARGET_VOLUME, FLOW_RATE)
start = time.perf_counter()

# Simple loop to capture volume and pressure data:
while pump.is_moving():
    volume_time.append(time.perf_counter() - start)
    volume_data.append(pump.get_infused_vol())
    force_time.append(time.perf_counter() - start)
    force_data.append(load_cell.get_force())
    pressure_time.append(time.perf_counter() - start)
    pressure_data.append(pressure_sensor.get_pressure())
    print(f"Pressure:\t {pressure_data[-1]} kPa")
    if pressure_data[-1] > MAXIMUM_PRESSURE:
        pump.stop()
        print(f"Pressure Exceeded {MAXIMUM_PRESSURE} kPa.")
        break
    # pump.set_flow_rate(new_calculated_flow_rate)
    time.sleep(0.01)

# Write data to output.csv in rows:
data = [volume_data, volume_time, pressure_data, pressure_time, force_time, force_data]
file = open(FILENAME, 'w', newline ='')
with file:
    write = csv.writer(file)
    write.writerows(data)
file.close()