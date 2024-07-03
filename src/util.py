"""Utility module."""

DEFAULT_PAD_LENGTH = 30


def pad_kv_pair_str(key: str, value, length: int = DEFAULT_PAD_LENGTH) -> str:
    """Created a padded key-value equivilency string. 

    Args:
        key (str): K-V key
        value: K-V value. 
        length (int): Total string length. Defaults to DEFAULT_PAD_LENGTH.

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
