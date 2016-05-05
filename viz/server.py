#!/usr/bin/env python

# Load required modules
import sys, os, json, argparse
import tornado.web, tornado.ioloop
from tornado.escape import json_encode

################################################################################
# Index page
class MainHandler(tornado.web.RequestHandler):
	def initialize(self, run_names):
		self.run_names = run_names

	def get(self):
		self.render("template.html", run_names=self.run_names)

# Data pages
class DataHandler(tornado.web.RequestHandler):
	def initialize(self, run_name, input_file):
		self.run_name = run_name
		with open(input_file, 'r') as IN:
			self.obj = IN.read()

	def get(self):
		self.write(self.obj)

################################################################################
# Argument parser
def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--input_files', type=str, nargs='*', required=True)
	parser.add_argument('-r', '--run_names', type=str, nargs='*', required=True)
	parser.add_argument('-p', '--port', type=int, required=False, default=8000)
	return parser

def run( args ):
	# Validate args
	assert( len(args.run_names) == len(args.input_files) )

	# Set up the server
	routes = [ (r"/", MainHandler, dict(run_names=args.run_names)),
	 		   (r'/bower_components/(.*)', tornado.web.StaticFileHandler, {'path': "bower_components"})]
	for run_name, input_file in zip(args.run_names, args.input_files):
		route = r"/data/{}.json".format(run_name)
		routes.append( (route, DataHandler, dict(run_name=run_name, input_file=input_file)) )

	# Start server
	app = tornado.web.Application(routes)
	app.listen(args.port)
	print 'Listening on port {}'.format(args.port)
	tornado.ioloop.IOLoop.current().start()

if __name__ == '__main__': run( get_parser().parse_args(sys.argv[1:]) )
