# Models

Talos uses [Pydantic](https://docs.pydantic.dev/latest/) models extensively, each time a file is read from, or written
to, disk. This validation means that we are always sure of object integrity.
