"""Index Page containing all routes for phactor."""

import os
import flask
from flask import request, session, redirect
import webkit


@webkit.app.route('/', methods=['GET', 'POST'])
def show_index():
    """Route for index '/'."""
    context = {}
    session.clear()
    return flask.render_template("index.html", **context)


@webkit.app.route('/favicon.ico')
def favicon():
    """Route for favicon."""
    return flask.send_from_directory(
                os.path.join(webkit.app.root_path, 'static'), 'favicon.ico')


@webkit.app.route('/invimg/<id>')
def show_invimg(id):
    """Route for inv image."""
    loc = webkit.app.config['INV_IMGS_FOLDER']

    return flask.send_from_directory(
                os.path.join(loc), str(id)+'.png')


