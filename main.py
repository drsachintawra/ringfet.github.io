from model import InputForm                              #//////////import modules
from flask import Flask, render_template, request
from demo3 import *
app = Flask(__name__)                                    # //////////create applicaton


@app.route('/', defaults={'page': 'index'})              #///////////////URL assign
@app.route('/<page>.html',  methods=['GET', 'POST'])
def index(page):

    form = InputForm(request.form)                          #data request and validation from WTforms module
    if request.method == 'POST' and form.validate():
            ra1,ra2,ra3,ra4,ra5,ra6 ,ra7,ra73,ra71,ra72,ra61,ra62,ra8,ra81,ra9,ra91 =  ash(form.Vgs.data,form.Vds.data) #,form.T.data,form.tox.data,form.Na.data,form.Nd.data'')  # return output from main file
    else:
        ra1 = None
        ra2 = None
        ra3= None
        ra4 = None
        ra5= None
        ra6 = None
        ra7 = None
        ra73 = None
        ra71 ,ra72 = None , None
        ra61 , ra62 = None , None
        ra8 , ra9 = None , None
        ra81 , ra91 = None , None

    return render_template('{}.html'.format(page), form=form,r1=ra1 , r2=ra2 ,r3=ra3,r4=ra4,r5=ra5 , r6 = ra6 ,
                           r7 = ra7 , r73 = ra73 ,r71 = ra71,r72 = ra72 , r61 = ra61,r62 = ra62 ,
                           r8 = ra8,r81 = ra81 , r9 = ra9, r91 = ra91)

if __name__ == '__main__':                      #///////////  calling the applicaton by name
    app.debug = True
    app.run(host = '0.0.0.0',port=5000)