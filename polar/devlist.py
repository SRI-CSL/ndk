#!/usr/bin/env python

import asyncio
from bleak import BleakScanner

all = []

async def find_devices():
    devices = await BleakScanner.discover()
    all = devices
    h10_list = []
    for d in devices:
        if d.name[0:9] == 'Polar H10':
            h10_list += [d]
            print(d.address, d.name, d.rssi)
    return h10_list

loop = asyncio.get_event_loop()
x = loop.run_until_complete(find_devices())

print(x)

# 293C33D4-2140-4D68-B85A-F5495CDA6D6C: Polar H10 78F5B92F
# 61819F2C-AE8F-4F05-BB0F-E9A541D103A7: Polar H10 78F5B92F 78:F5:B9:2F
# 61819F2C-AE8F-4F05-BB0F-E9A541D103A7: Polar H10 78F5B92F
