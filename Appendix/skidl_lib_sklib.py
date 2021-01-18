from skidl import Pin, Part, Alias, SchLib, SKIDL, TEMPLATE

SKIDL_lib_version = '0.0.1'

skidl_lib = SchLib(tool=SKIDL).add_parts(*[
        Part(**{ 'name':'EXPI', 'dest':TEMPLATE, 'tool':SKIDL, '_match_pin_regex':False, 'rise_delay_time':5, '_aliases':Alias({'exponentialcurrent', 'expi', 'EXPONENTIALCURRENT'}), 'fall_delay_time':5, 'keywords':'exponential current source', 'fall_time_constant':5, 'rise_time_constant':5, 'description':'Exponential current source', 'pyspice':{'name': 'ExponentialCurrentSource', 'kw': {'initial_value': 'initial_value', 'pulsed_value': 'pulsed_value', 'rise_delay_time': 'rise_delay_time', 'rise_time_constant': 'rise_time_constant', 'fall_delay_time': 'fall_delay_time', 'fall_time_constant': 'fall_time_constant', 'p': 'node_plus', 'n': 'node_minus'}, 'add': <function add_part_to_circuit at 0x7f42f18f6170>}, 'initial_value':5, 'pulsed_value':5, 'ref_prefix':'I', 'num_units':1, 'fplist':None, 'do_erc':True, 'aliases':Alias({'exponentialcurrent', 'expi', 'EXPONENTIALCURRENT'}), 'pin':None, 'footprint':None, 'pins':[
            Pin(num='1',name='p',func=Pin.types.PASSIVE,do_erc=True),
            Pin(num='2',name='n',func=Pin.types.PASSIVE,do_erc=True)] })])