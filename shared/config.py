"""
Config
"""

import os


APPLICATION_ROOT = '/'

SECRET_KEY = b'xe1\x14:\x92i\x99:\xb2:U\xfb\xa9h\xc8<\xe0\xd3\x17\xba\xf8R6\xfd\xd9'  # noqa: E501  pylint: disable=line-too-long
SESSION_COOKIE_NAME = 'login'


DATA_FOLDER = os.path.join(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
    'shared/data',
)

MODELS_FOLDER = os.path.join(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__))),
    'shared/models',
)




ALLOWED_EXTENSIONS = set(['png', 'jpg', 'jpeg', 'gif'])

MAX_CONTENT_LENGTH = 16 * 1024 * 1024

SEND_FILE_MAX_AGE_DEFAULT = 0
