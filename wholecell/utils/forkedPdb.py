import sys

try:
	from IPython.core.debugger import Pdb as _IpyPdb

	class ForkedPdb(_IpyPdb):
		""" Usage:
			==========
			from wholecell.utils import ForkedPdb
			ForkedPdb().set_trace()
			==========
			Don't forget to actually instantiate an
			instance of ForkedPdb when calling set_trace()

			Note: Python 3 has a global function breakpoint().
		"""
		def set_trace(self, frame = None):
			_stdin = sys.stdin
			sys.stdin = open("/dev/stdin")
			if frame is None:
				# noinspection PyUnresolvedReferences
				frame = sys._getframe().f_back
			_IpyPdb("Linux").set_trace(frame)

except ImportError:
	class ForkedPdb:  # type: ignore[no-redef]
		"""Stub — IPython is not installed.  Install it to use ForkedPdb."""
		def set_trace(self, frame = None):
			raise ImportError(
				"ForkedPdb requires IPython. Install it with: pip install ipython"
			)
