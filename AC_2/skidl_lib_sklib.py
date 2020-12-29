from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'SINEV', 'dest':TEMPLATE, 'tool':SKIDL, 'dc_value':UnitValue(10 V), '_aliases':Alias({'SINUSOIDALVOLTAGE', 'sinev', 'sinusoidalvoltage'}), '_match_pin_regex':False, 'keywords':'sinusoidal voltage source', 'ac_magnitude':UnitValue(10 V), 'pyspice':{'name': 'SinusoidalVoltageSource', 'kw': {'dc_offset': 'dc_offset', 'ac_magnitude': 'ac_magnitude', 'ac_phase': 'ac_phase', 'offset': 'offset', 'amplitude': 'amplitude', 'frequency': 'frequency', 'delay': 'delay', 'damping_factor': 'damping_factor', 'p': 'node_plus', 'n': 'node_minus'}, 'add': <function add_part_to_circuit at 0x7f6974541ef0>}, 'description':'Sinusoidal voltage source', 'ref_prefix':'V', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'SINUSOIDALVOLTAGE', 'sinev', 'sinusoidalvoltage'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'C', 'dest':TEMPLATE, 'tool':SKIDL, '_aliases':Alias({'CAP', 'cap'}), '_match_pin_regex':False, 'keywords':'cap capacitor', 'pyspice':{'name': 'C', 'kw': {'value': 'capacitance', 'capacitance': 'capacitance', 'model': 'model', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'ic': 'initial_condition', 'initial_condition': 'initial_condition', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f6974541ef0>}, 'description':'Capacitor', 'ref_prefix':'C', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'CAP', 'cap'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'R', 'dest':TEMPLATE, 'tool':SKIDL, 'keywords':'res resistor', 'description':'Resistor', '_match_pin_regex':False, 'pyspice':{'name': 'R', 'kw': {'value': 'resistance', 'resistance': 'resistance', 'ac': 'ac', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'noisy': 'noisy', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f6974541ef0>}, 'ref_prefix':'R', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias(), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])