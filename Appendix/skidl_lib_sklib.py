from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'C', 'dest':TEMPLATE, 'tool':SKIDL, 'dtemp':5, 'm':5, 'pyspice':{'name': 'C', 'kw': {'value': 'capacitance', 'capacitance': 'capacitance', 'model': 'model', 'multiplier': 'multiplier', 'm': 'multiplier', 'scale': 'scale', 'temp': 'temperature', 'temperature': 'temperature', 'dtemp': 'device_temperature', 'device_temperature': 'device_temperature', 'ic': 'initial_condition', 'initial_condition': 'initial_condition', 'p': 'plus', 'n': 'minus'}, 'add': <function add_part_to_circuit at 0x7efd40622dd0>}, 'description':'Capacitor', 'keywords':'cap capacitor', 'temp':5, 'scale':5, 'ic':5, '_match_pin_regex':False, '_aliases':Alias({'CAP', 'cap'}), 'ref_prefix':'C', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'CAP', 'cap'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])