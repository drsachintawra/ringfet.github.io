from model import InputForm
from flask import Flask, render_template, request
from tcad import *
app = Flask(__name__)


@app.route('/', defaults={'page': 'index'})              #///////////////URL assign
@app.route('/<page>.html',  methods=['GET', 'POST'])
def index(page):

    form = InputForm(request.form)
    if request.method == 'POST' and form.validate():
            ra1,ra2 =  as0()
    else:
        ra1 = None
        ra2 = None


    return render_template('{}.html'.format(page), form=form,r1=ra1 , r2=ra2 )


if __name__ == '__main__':
    app.run(debug=True)