from docutils import nodes
from docutils.parsers.rst import directives

CODE = """

<p><b><a href="javascript:doNothing()" onClick="popupHandle=popup(%(name)s)"><img src="%(depth)s" width=60 height=60 align=left>Watch the Video</a></b></p>
</br>

"""

PARAM = """\n    <param name="%s" value="%s"></param>"""

def popup(name, args, options, content, lineno,
            contentOffset, blockText, state, stateMachine):
    """ Restructured text extension for popups """
    if len(content) == 0:
        return
    string_vars = {
        'name': content[0],
        'depth': content[1],
        'extra': ''
        }
    #import ipdb; ipdb.set_trace()
    extra_args = content[1:] # Because self.content[0] is ID
    extra_args = [ea.strip().split("=") for ea in extra_args] # key=value
    extra_args = [ea for ea in extra_args if len(ea) == 2] # drop bad lines
    extra_args = dict(extra_args)
    if 'depth' in extra_args:
        string_vars['depth'] = extra_args.pop('depth')
    if extra_args:
        params = [PARAM % (key, extra_args[key]) for key in extra_args]
        string_vars['extra'] = "".join(params)
    return [nodes.raw('', CODE % (string_vars), format='html')]
popup.content = True
directives.register_directive('popup', popup)


def setup(foo):
    pass

