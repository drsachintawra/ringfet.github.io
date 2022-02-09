from wtforms import Form, FloatField, validators


class InputForm(Form):
    Vgs = FloatField(
        label='(V)', default=0.0,
        validators=[validators.InputRequired()])
    Vds = FloatField(
        label='(V)', default=0.0,
        validators=[validators.InputRequired()])
'''    T = FloatField(
        label='(K)', default=300,
        validators=[validators.InputRequired()])
    tox = FloatField(
        label='(m)', default=2e-9,
        validators=[validators.InputRequired()])
    Na = FloatField(
        label='(cm-3)', default=1e26,
        validators=[validators.InputRequired()])
    Nd = FloatField(
        label='(cm-3)', default=2e18,
        validators=[validators.InputRequired()])'''
