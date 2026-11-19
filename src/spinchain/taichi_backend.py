"""Optional Vulkan probe. Production remains on validated float64 Numba unless parity passes."""
def probe_vulkan():
    try:
        import taichi as ti
        return {"installed":True,"version":str(ti.__version__),"vulkan_compiled":bool(ti._lib.core.with_vulkan()),"f64_required":True,"approved_for_production":False}
    except Exception as e:return {"installed":False,"error":str(e),"approved_for_production":False}
