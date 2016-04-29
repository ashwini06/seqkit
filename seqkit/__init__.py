"""Make this a module"""
import yaml

# version number of the package
__version__ = "1.0"

# config to be used later
CONFIG = {}

def load_yaml_config(config_file):
    """Load YAML config file
    :param str config_file: The path to the configuration file.
    :returns: A dict of the parsed config file.
    :rtype: dict
    :raises IOError: If the config file cannot be opened.
    """
    if type(config_file) is file:
        CONFIG.update(yaml.load(config_file) or {})
        return CONFIG
    else:
        try:
            with open(config_file, 'r') as f:
                CONFIG.update(yaml.load(f))
                return CONFIG
        except IOError as e:
            e.message = "Could not open configuration file \"{}\".".format(config_file)
            raise e