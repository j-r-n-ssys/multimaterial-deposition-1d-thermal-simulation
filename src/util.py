from typing import Any


def pad_kv_pair_str(key: str, value, length: int = 30) -> str:
    """print a padded 

    Args:
        key (str): _description_
        value (_type_): _description_
        length (int): Total string length.

    Returns:
        str: _description_
    """

    if not isinstance(length, int):
        raise TypeError('Pad length must be an integer.')
    elif length <= 0:
        raise ValueError('Pad length must be greater than zero.')

    remainder = (length - 1) - len(key)

    if remainder <= 0:
        return f'{key} = {value}'
    else:
        return f'{key}{' '*remainder} = {value}'
