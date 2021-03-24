#!/usr/bin/env python

import os
import sys
import math
import time
import wfdb
import shutil
import signal
import asyncio
import os.path
import numpy as np
import nde.gen
import nde.ssread
import nde.stim as stim
import simpleaudio
import scipy.io.wavfile
from ndk.ds import wfdb_to_nbf

import pandas as pd
from bleak import BleakClient
from bleak.uuids import uuid16_dict
import matplotlib.pyplot as plt
import matplotlib
import simpleaudio as sa
import sounddevice as sd

matplotlib.use('TkAgg')
PLOT_IT = False

if (len(sys.argv) > 1):
    record_name = sys.argv[1]
    print("Sending data to {}".format(sys.argv[1]))

if (len(sys.argv) > 2):
    if sys.argv[2].isdigit():
        nminutes = int(sys.argv[2])
        nsecs = nminutes * 60
        print("Recording for {} minutes (= {} seconds = {} ticks).".format(nminutes, nsecs, nsecs*130))
        wave_buffer = None
    else:
        ssfile = sys.argv[2]
        state_machine = nde.ssread.ssread(ssfile)
        nsecs = state_machine.duration
        nminutes = nsecs / 60
        print("Recording for {} minutes (= {} seconds = {} ticks).".format(nminutes, nsecs, nsecs*130))
        seq = nde.gen.gen_sequence(state_machine)
        wavfile = ssfile[:-3] + ".wav"
        print("Generating wav file for stimulus in {}".format(wavfile))
        stream = stim.WAV_out(wavfile)
        tprev = 0
        sampnum = 0
        for i in range(len(seq)):
            time, amp, wave = seq[i]
            if i > 0:
                delta = int ( stream.samprate * (time - tprev) )
            else:
                delta = 0
            sampnum = stim.put_spike( stream, wave, delta, sampnum, scale=amp)
            tprev = time
   
        stream.put_sample(0.0, finish=True)
        stream.close()
        print("...Done.")
        sr,wave_buffer = scipy.io.wavfile.read(wavfile)
        wave_buffer = wave_buffer.reshape( ( len(wave_buffer), 1) )
#        waveobj = sa.WaveObject.from_wave_file(wavfile)

""" Predefined UUID (Universal Unique Identifier) mapping are based on Heart Rate GATT service Protocol that most
Fitness/Heart Rate device manufacturer follow (Polar H10 in this case) to obtain a specific response input from 
the device acting as an API """

uuid16_dict = {v: k for k, v in uuid16_dict.items()}

## This is the device MAC ID, please update with your device ID
# ADDRESS = "D4:52:48:88:EA:04"
ADDRESS = "61819F2C-AE8F-4F05-BB0F-E9A541D103A7"

## UUID for model number ##
MODEL_NBR_UUID = "0000{0:x}-0000-1000-8000-00805f9b34fb".format(
    uuid16_dict.get("Model Number String")
)


## UUID for manufacturer name ##
MANUFACTURER_NAME_UUID = "0000{0:x}-0000-1000-8000-00805f9b34fb".format(
    uuid16_dict.get("Manufacturer Name String")
)

## UUID for battery level ##
BATTERY_LEVEL_UUID = "0000{0:x}-0000-1000-8000-00805f9b34fb".format(
    uuid16_dict.get("Battery Level")
)

## UUID for connection establsihment with device ##
PMD_SERVICE = "FB005C80-02E7-F387-1CAD-8ACD2D8DF0C8"

## UUID for Request of stream settings ##
PMD_CONTROL = "FB005C81-02E7-F387-1CAD-8ACD2D8DF0C8"

## UUID for Request of start stream ##
PMD_DATA = "FB005C82-02E7-F387-1CAD-8ACD2D8DF0C8"

## UUID for Request of ECG Stream ##
ECG_WRITE = bytearray([0x02, 0x00, 0x00, 0x01, 0x82, 0x00, 0x01, 0x01, 0x0E, 0x00])

## For Plolar H10  sampling frequency ##
ECG_SAMPLING_FREQ = 130

time_0 = -1
psec = 0
ecg_session_data = []
ecg_session_time = []


def save_signal(data, name):
    p_signal = np.asarray(data, dtype=np.float64)
    n = len(data)
    p_signal = p_signal.reshape((n,1))
    print(p_signal)
    record_name = name
    dat_file_name = ['{}.dat'.format(name)]
    units = ['qV']
    fmt = ['16']
    adc_gain = [1]
    baseline = [0]
    samps_per_frame = [1]

    rec = wfdb.io.wrsamp(name,
                         fs=130,
                         p_signal = p_signal,
                         units=units,
                         sig_name=['V1'],
                         fmt = fmt)
    return(rec)

def save_and_process(data, name, channels=['V1']):
    save_signal(data, name)
    src_hdr = name + '.hea'
    src_dat = name + '.dat'
    to_dir = './'+name+'/'
    print("Saving wfdb record {} to directory {}.".format(name, to_dir))
    wfdb_to_nbf(name, to_dir, ['V1'])
    os.system("cp {} {}".format(ssfile, to_dir))
    # shutil.copyfile(ssfile, to_dir)
    print("Saving a copy of wfdb header file {} to {}".format(src_hdr, to_dir))
    os.system("cp {} {}".format(src_hdr, to_dir))
    # shutil.copyfile(src_hdr, to_dir)
    print("Saving a copy of wfdb data file {} to {}".format(src_dat, to_dir))
    os.system("cp {} {}".format(src_dat, to_dir))
    # shutil.copyfile(src_dat, to_dir)
    print("Saving state machine events to file {}.".format(to_dir+'events.dat'))
    state_machine.save_events( to_dir+'events.dat' )


## Positoning/Pinnning the real-time plot window on the screen
def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    backend = matplotlib.get_backend()
    if backend == "TkAgg":
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == "WXAgg":
        f.canvas.manager.window.SetPosition((x, y))
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)


## Keyboard Interrupt Handler
def keyboardInterrupt_handler(signum, frame):
    global playback_event   # this sucks
    print("  key board interrupt received...")
    print("----------------Recording stopped------------------------")
    if playback_event:
        playback_event.set()
    else:
        # This is just a timed collect with no playback, so save and exit:
        save_and_process(ecg_session_data, record_name)
        # sys.exit(0)
        quit()



## Bit conversion of the Hexadecimal stream
def data_conv(sender, data):
    global psec, time_0
    if data[0] == 0x00:
        timestamp = convert_to_unsigned_long(data, 1, 8)
        step = 3
        samples = data[10:]
        offset = 0
        sec = int(timestamp / 1000000000)
        if time_0 == -1:
            time_0 = sec
        if sec > psec:
            print(sec-time_0, end=' ')
            sys.stdout.flush()
            psec = sec
        while offset < len(samples):
            ecg = convert_array_to_signed_int(samples, offset, step)
            offset += step
            ecg_session_data.extend([ecg])
            ecg_session_time.extend([timestamp])


def convert_array_to_signed_int(data, offset, length):
    return int.from_bytes(
        bytearray(data[offset : offset + length]), byteorder="little", signed=True,
    )


def convert_to_unsigned_long(data, offset, length):
    return int.from_bytes(
        bytearray(data[offset : offset + length]), byteorder="little", signed=False,
    )


 
async def play_buffer(buffer, **kwargs):
    if buffer is None:
        return None

    print("Setting up audio stream")
    loop = asyncio.get_event_loop()
    event = asyncio.Event()
    idx = 0

    def callback(outdata, frame_count, time_info, status):
        nonlocal idx
        if status:
            print(status)
        remainder = len(buffer) - idx
        if remainder == 0:
            loop.call_soon_threadsafe(event.set)
            raise sd.CallbackStop
        valid_frames = frame_count if remainder >= frame_count else remainder
        outdata[:valid_frames] = buffer[idx:idx + valid_frames]
        outdata[valid_frames:] = 0
        idx += valid_frames

    print("Callback defined.  Opening audio stream (type={})".format(buffer.dtype))
    stream = sd.OutputStream(samplerate=sr, callback=callback, dtype=buffer.dtype,
                             channels=1, blocksize=8192, **kwargs)
    print("Stream opened: ", stream)
    
    # We will want to await this to wait for playback to reach the end of the buffer:
    return event, stream
# We want this to run concurrently
#    with stream:
#        await event.wait()

def play_buffer2(buffer, **kwargs):
    if buffer is None:
        return None

    print("Setting up audio stream")
    idx = 0

    def callback(outdata, frame_count, time_info, status):
        nonlocal idx
        if status:
            print(status)
        remainder = len(buffer) - idx
        if remainder == 0:
            raise sd.CallbackStop()
        valid_frames = frame_count if remainder >= frame_count else remainder
        outdata[:valid_frames] = buffer[idx:idx + valid_frames]
        outdata[valid_frames:] = 0
        idx += valid_frames

    print("Callback defined.  Opening audio stream")
    stream = sd.OutputStream(samplerate=sr, callback=callback, dtype=buffer.dtype,
                             channels=1, blocksize=8192, **kwargs)
    print("Stream opened: ", stream)
    return stream



## Aynchronous task to start the data stream for ECG ##
async def run(client, debug=False):
    global playback_event
    ## Writing chracterstic description to control point for request of UUID (defined above) ##

    # Should this look anything like the ECG stream?
    
    await client.is_connected()
    print("---------Device connected--------------")

    model_number = await client.read_gatt_char(MODEL_NBR_UUID)
    print("Model Number: {0}".format("".join(map(chr, model_number))))

    manufacturer_name = await client.read_gatt_char(MANUFACTURER_NAME_UUID)
    print("Manufacturer Name: {0}".format("".join(map(chr, manufacturer_name))))

    battery_level = await client.read_gatt_char(BATTERY_LEVEL_UUID)
    print("Battery Level: {0}%".format(int(battery_level[0])))

    att_read = await client.read_gatt_char(PMD_CONTROL)

    await client.write_gatt_char(PMD_CONTROL, ECG_WRITE)

    ## ECG stream started
    await client.start_notify(PMD_DATA, data_conv)

    print("Collecting ECG data...")

    print("Starting stimulation playback")
    playback_event, playback_stream = await play_buffer(wave_buffer)

    n = ECG_SAMPLING_FREQ

    ## Collecting ECG data for nsecs seconds
    if playback_event is None:
        await asyncio.sleep(nsecs)
        print("...")
        print("We ran for {} seconds.  Stopping collection.".format(nsecs))
    else:
        # End of playback == end of experiment:
        with playback_stream:
            await playback_event.wait()
        playback_stream.close()
        print("...")
        print("Reached the end of the stimulus program.  Stopping collection.")

    ## Stop the stream once data is collected
    await client.stop_notify(PMD_DATA)
    print("Stopping ECG data...")
    print("[CLOSED] application closed.")

    save_and_process(ecg_session_data, record_name)

    quit()
    # sys.exit(0)


async def main():
    try:
        async with BleakClient(ADDRESS) as client:
            signal.signal(signal.SIGINT, keyboardInterrupt_handler)
            tasks = [
                asyncio.ensure_future(run(client, True)),
                
            ]

            await asyncio.gather(*tasks)
    except:
        pass


if __name__ == "__main__":
    os.environ["PYTHONASYNCIODEBUG"] = str(1)
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)
    loop.run_until_complete(main())
