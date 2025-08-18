from logging import log
import logging.config
import atexit


def setup_logging(args) -> None:

    logging_level = "DEBUG" if args.debug else "INFO"
    
    config = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "simple": {
                "format": "%(message)s"
            }
        },
        "handlers": {
            "stdout": {
                "class": "logging.StreamHandler",
                "formatter": "simple",
                "stream": "ext://sys.stdout"
            },
            "file": {
                "class": "logging.FileHandler",
                "formatter": "simple",
                "filename": args.output_file
            },
            "queue_handler": {
                "class": "logging.handlers.QueueHandler",
                "handlers": [ "stdout", "file" ],
                "respect_handler_level": True
            }
        },
        "loggers": {
            "root": {"level": logging_level, "handlers": ["queue_handler"]},
            "script": {"level": logging_level}
        }
    }

    logging.config.dictConfig(config)
    queue_handler = logging.getHandlerByName("queue_handler")
    if queue_handler is not None:
        queue_handler.listener.start()
        atexit.register(queue_handler.listener.stop)

    logging.getLogger("matplotlib").setLevel(logging.CRITICAL)

