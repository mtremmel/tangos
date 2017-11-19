from __future__ import absolute_import
from tangos.core import extraction_patterns
from . import BuiltinFunction
from .. import StoredProperty, FixedInput, LiveProperty


@BuiltinFunction.register
def raw(halos, values):
    return values
#raw.set_input_options(0, assert_class=StoredProperty)

@raw.set_initialisation
def raw_initialisation(input):
    print "in raw init", isinstance(input,LiveProperty), input
    if isinstance(input, LiveProperty):
        print "setting eval pattern"
        input.set_evaluation_pattern('_evaluate_function')
    else:
        input.set_extraction_pattern(extraction_patterns.halo_property_raw_value_getter)


@BuiltinFunction.register
def reassemble(halos, values, *options):
    return values
#reassemble.set_input_options(0, assert_class=StoredProperty)

@reassemble.set_initialisation
def reassemble_initialisation(input, *options):
    print "in reassemble init", isinstance(input,LiveProperty), input, type(input)
    options_values = []
    for option in options:
        if isinstance(option, FixedInput):
            options_values.append(option.proxy_value())
        else:
            raise TypeError("Options to 'reassemble' must be fixed numbers or strings")
    if isinstance(input, LiveProperty):
        print "setting eval pattern"
        input.set_evaluation_pattern('_evaluate_function_with_reassemble', *options_values)
    else:
        input.set_extraction_pattern(extraction_patterns.HaloPropertyValueWithReassemblyOptionsGetter(*options_values))
