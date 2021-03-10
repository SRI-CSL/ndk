#!/usr/bin/env python

import asyncio
from bleak import BleakScanner

all = []

async def run():
    devices = await BleakScanner.discover()
    all = devices
    for d in devices:
        print(d)

loop = asyncio.get_event_loop()
loop.run_until_complete(run())


# 61819F2C-AE8F-4F05-BB0F-E9A541D103A7: Polar H10 78F5B92F 78:F5:B9:2F
# 61819F2C-AE8F-4F05-BB0F-E9A541D103A7: Polar H10 78F5B92F
