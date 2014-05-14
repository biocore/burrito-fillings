#!/usr/bin/env python

"""
This module provides pass-through access to PyCogent's denoiser code. It's a
bit of a hack, but it allows us to remove the direct dependency on PyCogent by
centralizing the denoiser code with all of the other PyCogent code that is
targeted either for complete re-write or removal pending benchmarks. The basic
idea is that it's not worth porting this code anywhere now because it's days are
numbered, but we still need to be able to access it for the time being.

May The Biological Research (application controller) Obsolescence Kit (brokit)
be as short-lived as its name is awesome.

"""

from cogent.parse.flowgram import (Flowgram, build_averaged_flowgram,
                                   seq_to_flow)
from cogent.parse.flowgram_parser import lazy_parse_sff_handle, get_header_info
from cogent.parse.flowgram_collection import (FlowgramCollection, parse_sff)
