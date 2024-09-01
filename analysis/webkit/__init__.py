"""
webkit package initializer.
"""

import flask

app = flask.Flask(__name__)  # pylint: disable=invalid-name

app.config.from_object("shared.config")
app.config.from_envvar('APP_SETTINGS', silent=True)

SESSION_TYPE = 'filesystem'
app.config.from_object(__name__)

import webkit.views  # noqa: E402  pylint: disable=wrong-import-position
import webkit.api  # noqa: E402  pylint: disable=wrong-import-position