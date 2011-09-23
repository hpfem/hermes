from docutils import nodes
from docutils.parsers.rst import directives

CODE = """\
<iframe class="youtube-player" type="text/html" width="640" height="385" src="http://www.youtube.com/embed/%(yid)s" frameborder="0">
</iframe>
"""

PARAM = """\n    <param name="%s" value="%s"></param>"""

def youtube(name, args, options, content, lineno,
            contentOffset, blockText, state, stateMachine):
    """ Restructured text extension for inserting youtube embedded videos """
    if len(content) == 0:
        return
    string_vars = {
        'yid': content[0],
        'width': 528,
        'height': 400,
        'extra': ''
        }
    extra_args = content[1:] # Because content[0] is ID
    extra_args = [ea.strip().split("=") for ea in extra_args] # key=value
    extra_args = [ea for ea in extra_args if len(ea) == 2] # drop bad lines
    extra_args = dict(extra_args)
    if 'width' in extra_args:
        string_vars['width'] = extra_args.pop('width')
    if 'height' in extra_args:
        string_vars['height'] = extra_args.pop('height')
    if extra_args:
        params = [PARAM % (key, extra_args[key]) for key in extra_args]
        string_vars['extra'] = "".join(params)
    return [nodes.raw('', CODE % (string_vars), format='html')]
youtube.content = True
directives.register_directive('youtube', youtube)


def setup(foo):
    pass
