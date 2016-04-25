"""Make this a module"""
import yaml

#version number of the package
__version__ = "1.0"

def load_yaml_config(config_file):
    """Load YAML config file
    :param str config_file: The path to the configuration file.
    :returns: A dict of the parsed config file.
    :rtype: dict
    :raises IOError: If the config file cannot be opened.
    """
    if type(config_file) is file:
        return yaml.load(config_file) or {}
    else:
        try:
            with open(config_file, 'r') as f:
                return yaml.load(f)
        except IOError as e:
            e.message = "Could not open configuration file \"{}\".".format(config_file)
            raise e