"""Data adapters to make available input data from different sources.

This module implements classes for the parsing of ortholog sequence data.

"""

class HamstrAdapter:
	"""Parse HaMStR output.
	"""
	dir = ''

	def __init__(self, dir):
		self.dir = dir

	def __repr__(self):
		return "<HamstrAdapter(dir='%s')>" % self.dir