# Frequently asked questions (FAQ)


### Which tractograms can I use with COMMIT?

Any tractogram can be fed to the COMMIT framework, as long as each tract is represented as a polyline, i.e. sequence of consecutive segments. At the moment, however, COMMIT only reads tractograms in the [.TRK](http://www.trackvis.org/docs/?subsect=fileformat) and [.TCK](https://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html#tracks-file-format-tck) file formats. For other formats, several converters are available, e.g. [DIPY](http://dipy.org).

### Where is the old MATLAB version?

The old MATLAB version is still available [as a tag in the repository](https://github.com/daducci/COMMIT/releases/tag/MATLAB). Please note, however, that this version is no longer mantained.

### Why this transition to Python?

We decided to re-implement our tool in Python with the aim to be more compatible with existing tools and libraries in the field, in particular [DIPY](http://dipy.org), and allow an easier integration in existing pipelines. 
