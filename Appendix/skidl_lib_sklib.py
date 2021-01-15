from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'C', 'dest':TEMPLATE, 'tool':SKIDL, 'keywords':'cap capacitor', 'pyspice':{'name': 'C', 'kw': {'value': 'capacitance', 'capacitance': 'capacitance', 'model': 'model', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'ic': 'initial_condition', 'initial_condition': 'initial_condition', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7f7e7803c170>}, 'm':5, 'temp':5, '_aliases':Alias({'cap', 'CAP'}), 'scale':5, 'ic':5, '_match_pin_regex':False, 'dtemp':5, 'description':'Capacitor', 'ref_prefix':'C', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'cap', 'CAP'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])