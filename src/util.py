"""Utility module."""


def pad_kv_pair_str(key: str, value, length: int = 30) -> str:
    """Created a padded key-value equivilency string. 

    Args:
        key (str): K-V key
        value: K-V value. 
        length (int): Total string length.

    Returns:
        str: Padded key-value string.
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
