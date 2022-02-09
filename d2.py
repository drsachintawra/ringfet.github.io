from flask import Flask
app = Flask(__name__)

@app.route('/')
def hello_world():
   return "Hey, I am Ashish"

if __name__ == '__main__':
   app.run()

'''from flask import Flask, render_template, request
app = Flask(__name__)



@app.route('/', defaults={'page': 'home'})


@app.route('/<page>.html',  methods=['GET', 'POST'])
def index(page):
    return render_template('{}.html'.format(page))


@app.route('/', defaults= {'page' : 'home'})
@app.route('/<page>.html',  methods=['GET', 'POST'])
def index(page):
    return render_template('{}.html'.format(page))

@app.route('/layout/<log>.html' , methods = ['GET' , 'POST' ])
def index1(log):
    return render_template('{}.html'.format(log))

if __name__ == '__main__':
    app.run(debug = True)
    #app.debug = True
    #app.run(host = '0.0.0.0',port=5000)/'''

