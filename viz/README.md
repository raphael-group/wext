# Visualization #

We provide scripts for viewing the output of the `compute_exclusivity.py` script in an interactive web application.

#### Setup ###

The visualization scripts have additional dependencies that must be installed.

1. **Python**. The web application is run using a [Tornado web server](http://www.tornadoweb.org/en/stable/). Install the required modules in `requirements.txt` using the following command:

        pip install -r requirements.txt

2. **Node, NPM, and Bower**. We use [Bower](https://bower.io) to manage Javascript and CSS dependencies. You can install Bower using the [Node Package Manager (NPM)](https://www.npmjs.com/), which requires [Node.js](https://nodejs.org/en/) to be installed. Once NPM is installed, install Bower and then use Bower to install the dependencies:

        npm install -g bower
        bower install

 If successful, this will create a directory `bower_components`.

### Usage ###

First, use `generate_viz_data.py` to generate a JSON file(s) which will serve as the input for the web server.

Then, run `server.py` to start the web server. By default, it will serve the web application port 8000, which you can view by pointing your browser to `http://localhost:8000`.


For additional details on usage and the arguments to the two scripts, please see [the wiki](https://github.com/raphael-group/weighted-exclusivity-test/wiki/VIsualization).
