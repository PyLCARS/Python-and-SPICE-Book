from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'SINEV', 'dest':TEMPLATE, 'tool':SKIDL, 'ac_magnitude':UnitValue(1 V), 'pyspice':{'name': 'SinusoidalVoltageSource', 'kw': {'dc_offset': 'dc_offset', 'ac_magnitude': 'ac_magnitude', 'ac_phase': 'ac_phase', 'offset': 'offset', 'amplitude': 'amplitude', 'frequency': 'frequency', 'delay': 'delay', 'damping_factor': 'damping_factor', 'p': 'node_plus', 'n': 'node_minus'}, 'add': <function add_part_to_circuit at 0x7f2850134170>}, 'keywords':'sinusoidal voltage source', '_match_pin_regex':False, 'description':'Sinusoidal voltage source', '_aliases':Alias({'SINUSOIDALVOLTAGE', 'sinusoidalvoltage', 'sinev'}), 'ref_prefix':'V', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'SINUSOIDALVOLTAGE', 'sinusoidalvoltage', 'sinev'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'R', 'dest':TEMPLATE, 'tool':SKIDL, 'description':'Resistor', 'pyspice':{'name': 'R', 'kw': {'value': 'resistance', 'resistance': 'resistance', 'ac': 'ac', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'noisy': 'noisy', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f2850134170>}, 'keywords':'res resistor', '_match_pin_regex':False, 'ref_prefix':'R', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias(), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'A', 'dest':TEMPLATE, 'tool':SKIDL, 'model':<skidl.tools.spice.XspiceModel object at 0x7f282d5dd390>, 'pyspice':{'name': 'A', 'kw': {'model': 'model'}, 'add': <function add_xspice_to_circuit at 0x7f285014df80>}, 'keywords':'XSPICE', '_match_pin_regex':False, 'description':'XSPICE code module', '_aliases':Alias({'xspice', 'XSPICE'}), 'ref_prefix':'A', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'xspice', 'XSPICE'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(name='io1',func=Pin.types.UNSPEC,do_erc=True),
            Pin(num=1,name='io2',func=Pin.types.UNSPEC,do_erc=True),
            Pin(num=2,name='mmf1',func=Pin.types.UNSPEC,do_erc=True),
            Pin(num=3,name='mmf2',func=Pin.types.UNSPEC,do_erc=True)] }),
        Part(**{ 'name':'V', 'dest':TEMPLATE, 'tool':SKIDL, 'pyspice':{'name': 'V', 'kw': {'value': 'dc_value', 'dc_value': 'dc_value', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f2850134170>}, 'keywords':'voltage source', '_match_pin_regex':False, 'description':'Voltage source', 'dc_value':0, '_aliases':Alias({'VS', 'v', 'ammeter', 'AMMETER', 'vs'}), 'ref_prefix':'V', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'VS', 'v', 'ammeter', 'AMMETER', 'vs'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])