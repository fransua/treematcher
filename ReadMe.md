[![Gitter Chat](http://img.shields.io/badge/chat-online-brightgreen.svg)](https://gitter.im/etetoolkit/treematcher)
[![Build Status](https://travis-ci.org/etetoolkit/treematcher.svg?branch=master)](https://travis-ci.org/etetoolkit/treematcher)

To test:

python treematcher.py

## To Do
- use module six for keeping compatible syntax for python 2/3
- figure out what to do with @.children in searches.
```
    Example:
    Pattern = """
    ( '"Primates" in @.children[0].lineage ' )' @.support >= 0.9 ';
    """
    Pattern fails once it runs into a leaf. A leaf will return an empty list for children.
    You can't take the lineage of an empty list.
    Perhaps we should print the error as a warning and return False in these cases.

```
