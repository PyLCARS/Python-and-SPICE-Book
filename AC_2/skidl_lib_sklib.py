from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'SINEV', 'dest':TEMPLATE, 'tool':SKIDL, 'description':'Sinusoidal voltage source', 'amplitude':UnitValue(10 V), 'frequency':FrequencyValue(10 kHz), 'keywords':'sinusoidal voltage source', '_match_pin_regex':False, '_aliases':Alias({'sinev', 'SINUSOIDALVOLTAGE', 'sinusoidalvoltage'}), 'pyspice':{'name': 'SinusoidalVoltageSource', 'kw': {'dc_offset': 'dc_offset', 'ac_magnitude': 'ac_magnitude', 'ac_phase': 'ac_phase', 'offset': 'offset', 'amplitude': 'amplitude', 'frequency': 'frequency', 'delay': 'delay', 'damping_factor': 'damping_factor', 'p': 'node_plus', 'n': 'node_minus'}, 'add': <function add_part_to_circuit at 0x7ff34c038dd0>}, 'ref_prefix':'V', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'sinev', 'SINUSOIDALVOLTAGE', 'sinusoidalvoltage'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'C', 'dest':TEMPLATE, 'tool':SKIDL, 'description':'Capacitor', 'keywords':'cap capacitor', '_match_pin_regex':False, '_aliases':Alias({'cap', 'CAP'}), 'pyspice':{'name': 'C', 'kw': {'value': 'capacitance', 'capacitance': 'capacitance', 'model': 'model', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'ic': 'initial_condition', 'initial_condition': 'initial_condition', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7ff34c038dd0>}, 'ref_prefix':'C', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'cap', 'CAP'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'R', 'dest':TEMPLATE, 'tool':SKIDL, 'description':'Resistor', 'keywords':'res resistor', 'pyspice':{'name': 'R', 'kw': {'value': 'resistance', 'resistance': 'resistance', 'ac': 'ac', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'noisy': 'noisy', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7ff34c038dd0>}, '_match_pin_regex':False, 'ref_prefix':'R', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias(), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])