from __future__ import annotations
import json,os,platform,subprocess,sys

def powershell(cmd):
    try:return subprocess.check_output(['powershell','-NoProfile','-Command',cmd],text=True,stderr=subprocess.STDOUT,timeout=30).strip()
    except Exception as e:return f'unavailable: {e}'
def probe():
    info={"platform":platform.platform(),"python":sys.version,"processor":platform.processor(),"cpu_count":os.cpu_count()}
    if os.name=='nt':
        info['cpu']=powershell("Get-CimInstance Win32_Processor | Select Name,NumberOfCores,NumberOfLogicalProcessors,MaxClockSpeed | ConvertTo-Json -Compress")
        info['gpus']=powershell("Get-CimInstance Win32_VideoController | Select Name,AdapterRAM,DriverVersion | ConvertTo-Json -Compress")
        info['memory']=powershell("Get-CimInstance Win32_ComputerSystem | Select TotalPhysicalMemory | ConvertTo-Json -Compress")
    try:
        import taichi as ti; info['taichi']=ti.__version__;info['vulkan_available']=bool(ti._lib.core.with_vulkan())
    except Exception as e:info['taichi_error']=str(e);info['vulkan_available']=False
    return info
