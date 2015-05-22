# burrito-fillings changelog

## Version 0.1.1 (2015-05-22)

* Updated handling of temporary files to make better use of python ``tempfile.gettempdir()`` for some of the most widely used burrito fillings ([#61](https://github.com/biocore/burrito-fillings/pull/61), [#64](https://github.com/biocore/burrito-fillings/pull/64)).
* Fixed bug where swarm wrapper would silently ignore ``swarm`` failures ([#67](https://github.com/biocore/burrito-fillings/pull/67), [biocore/qiime#2014](https://github.com/biocore/qiime/issues/2014)).
* Added ``__version__`` to ``bfillings/__init__.py`` so that other python packages have access to the version number ([#54](https://github.com/biocore/burrito-fillings/issues/54)).

## Version 0.1.0 (2014-11-12)

Initial release.
