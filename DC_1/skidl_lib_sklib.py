from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'V', 'dest':TEMPLATE, 'tool':SKIDL, 'keywords':'voltage source', 'pyspice':{'name': 'V', 'kw': {'value': 'dc_value', 'dc_value': 'dc_value', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f4dafdcf170>}, 'description':'Voltage source', 'dc_value':UnitValue(10 V), '_aliases':Alias({'ammeter', 'AMMETER', 'vs', 'v', 'VS'}), '_match_pin_regex':False, 'ref_prefix':'V', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'ammeter', 'AMMETER', 'vs', 'v', 'VS'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] }),
        Part(**{ 'name':'R', 'dest':TEMPLATE, 'tool':SKIDL, '_match_pin_regex':False, 'pyspice':{'name': 'R', 'kw': {'value': 'resistance', 'resistance': 'resistance', 'ac': 'ac', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'noisy': 'noisy', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f4dafdcf170>}, 'keywords':'res resistor', 'description':'Resistor', 'ref_prefix':'R', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias(), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])