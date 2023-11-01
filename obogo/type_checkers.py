from inspect import signature
from typing import Literal, get_args, get_origin
from decorator import decorator
import typing

@decorator
def literal_arg_checker(f, *args, **kwargs):
    sig        = signature(f)
    params     = sig.parameters
    bound_args = sig.bind(*args, **kwargs)
    
    for k_param, v_param in params.items():
        annot_type = get_origin(v_param.annotation)
        if annot_type is typing.Literal:
            allow_word = get_args(v_param.annotation)
            if k_param in bound_args.arguments:               
                bound_value = bound_args.arguments[k_param]
                if not  bound_value in allow_word:
                    raise TypeError(f"parameter {k_param} value \"{bound_value}\" is not a valid Literal {allow_word}")
    return f(*args, **kwargs)