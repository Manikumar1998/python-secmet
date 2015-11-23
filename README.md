SecMet - A Python Library to Handle Secondary Metabolite Annotations
====================================================================

This library is intended to make dealing with secondary metabolite annotations
as provided by databases like [MIBiG][mibig] or tools like
[antiSMASH][antismash] easier from within python programs.

At the moment, it just supports the GenBank file format used by [MIBiG][mibig]
and [antiSMASH][antismash], not the MIBiG JSON format. Support for this is
planned.


[mibig]: http://mibig.secondarymetabolites.org/
  "Minimum Information about a Biosynthetic Gene cluster (MIBiG) specification"

[antismash]: http://antismash.secondarymetabolites.org/
  "antibiotics & Secondary Metabolite Analysis Shell"


Installing
----------

SecMet is a normal Python module, so you can install it from a local repository
checkout using `python setup.py install`, or from PyPi using `pip install
secmet`.

Note that secmet is still in early development, and not fit for production use.


Packaging
---------

Because I perfer to work with markdown and Python forces me to use
reStructuredText, a Makefile helps in order to create a new PyPi release. The
`pandoc` utility is required. Run the Makefile using `make release` to create a
new secmet library sdist.

License
-------
    Copyright 2015  Kai Blin <kblin@biosustain.dtu.dk>

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
