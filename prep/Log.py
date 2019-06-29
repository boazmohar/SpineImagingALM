import functools
import inspect
import logging
import sys


def add_logging_to_file(logger_name, path, mode='a', level=logging.DEBUG,
                        log_format="%(name)s @ %(asctime)s - [%(levelname)s] %(module)s::%(funcName)s: %(message)s"):
    """
    
    :param logger_name: 
    :param path: 
    :param mode: 
    :param level: 
    :param log_format: 
    :return: 
    """
    logger = logging.getLogger(logger_name)
    fh = logging.handlers.RotatingFileHandler(path, maxBytes=(1048576 * 5), backupCount=7, mode=mode)
    fh.setFormatter(logging.Formatter(log_format))
    fh.setLevel(level)
    logger.addHandler(fh)


def name(item):
    " Return an item's name. "
    return item.__name__


def is_classmethod(instancemethod):
    " Determine if an instancemethod is a classmethod. "
    return instancemethod.im_self is not None


def is_class_private_name(name):
    " Determine if a name is a class private name. "
    # Exclude system defined names such as __init__, __add__ etc
    return name.startswith("__") and not name.endswith("__")


def method_name(method):
    """ Return a method's name.

    This function returns the name the method is accessed by from
    outside the class (i.e. it prefixes "private" methods appropriately).
    """
    mname = name(method)
    if is_class_private_name(mname):
        mname = "_%s%s" % (name(method.im_class), mname)
    return mname


def format_arg_value(arg_val, threshold=1000, replace_with='big_value'):
    """ Return a string representing a (name, value) pair. If value is larger then threshold it is changed to
     the value of replace_with

    >>> format_arg_value(('x', (1, 2, 3)))
    'x=(1, 2, 3)'
    >>> format_arg_value(('a'), range(1000), threshold=100, replace_with='big_value')
    'a=big_value'
    """
    arg, val = arg_val
    if sys.getsizeof(val, threshold) >= threshold:
        val='big_value'
    return "%s=%r" % (arg, val)


def echo(fn, logger):
    """ Echo calls to a function.

    Returns a decorated version of the input function which "echoes" calls
    made to it by writing out the function's name and the arguments it was
    called with.
    """
    # Unpack function's arg count, arg names, arg defaults
    try:
        code = fn.func_code
    except AttributeError:
        code = fn.__code__
    argcount = code.co_argcount
    argnames = code.co_varnames[:argcount]
    try:
        fn_defaults = fn.func_defaults or list()
    except AttributeError:
        fn_defaults = fn.__defaults__ or list()
    argdefs = dict(zip(argnames[-len(fn_defaults):], fn_defaults))
    
    @functools.wraps(fn)
    def wrapped(*v, **k):
        # Collect function arguments by chaining together positional,
        # defaulted, extra positional and keyword arguments.
        positional = map(format_arg_value, zip(argnames, v))
        defaulted = [format_arg_value((a, argdefs[a]))
                     for a in argnames[len(v):] if a not in k]
        keyword = map(format_arg_value, k.items())
        nameless = v[argcount:]
        if nameless is not None:
            n = len(nameless)
            names = []
            for i in range(n):
                names += 'pos' + str(i)
            nameless = map(format_arg_value, zip(names, nameless))
        args = positional + defaulted + nameless + keyword
        logger.debug("%s(%s)\n" % (name(fn), ", ".join(args)))
        return fn(*v, **k)

    return wrapped


def echo_instancemethod(klass, method, logger):
    """ Change an instancemethod so that calls to it are echoed.

    Replacing a classmethod is a little more tricky.
    See: http://www.python.org/doc/current/ref/types.html
    """
    mname = method_name(method)
    never_echo = "__str__", "__repr__",  # Avoid recursion printing method calls
    if mname in never_echo:
        pass
    elif is_classmethod(method):
        setattr(klass, mname, classmethod(echo(method.im_func, logger)))
    else:
        setattr(klass, mname, echo(method, logger))


def echo_class(klass, logger):
    """ Echo calls to class methods and static functions
    """
    for _, method in inspect.getmembers(klass, inspect.ismethod):
        echo_instancemethod(klass, method, logger)
    for _, fn in inspect.getmembers(klass, inspect.isfunction):
        setattr(klass, name(fn), staticmethod(echo(fn, logger)))


def echo_module(mod, logger):
    """ Echo calls to functions and methods in a module.
    """
    try:
        functions = [x for x in inspect.getmembers(mod, inspect.isfunction) if
                     x[1].func_globals['__name__'] == mod.__name__]
    except AttributeError:
        functions = [x for x in inspect.getmembers(mod, inspect.isfunction) if
                     x[1].__globals__['__name__'] == mod.__name__]
    classes = [x for x in inspect.getmembers(mod, inspect.isclass) if
               x[1].__module__ == mod.__name__]
    for fname, fn in functions:
        setattr(mod, fname, echo(fn, logger))
    for _, klass in classes:
        echo_class(klass, logger)
