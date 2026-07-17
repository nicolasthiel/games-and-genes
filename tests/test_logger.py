import logging
import tempfile
import unittest
from pathlib import Path

from games_and_genes.logger import setup_logging


class TestLoggingSetup(unittest.TestCase):

    def setUp(self):
        self.root_logger = logging.getLogger()
        self.original_handlers = self.root_logger.handlers[:]
        self.original_level = self.root_logger.level
        for handler in self.root_logger.handlers[:]:
            self.root_logger.removeHandler(handler)
            handler.close()

    def tearDown(self):
        for handler in self.root_logger.handlers[:]:
            self.root_logger.removeHandler(handler)
            handler.close()
        for handler in self.original_handlers:
            self.root_logger.addHandler(handler)
        self.root_logger.setLevel(self.original_level)

    def test_setup_logging_is_idempotent(self):
        temp_dir = tempfile.TemporaryDirectory()
        try:
            logger = setup_logging({"level": "DEBUG", "log_to_file": True}, temp_dir.name)
            setup_logging({"level": "DEBUG", "log_to_file": True}, temp_dir.name)

            self.assertIs(logger, self.root_logger)
            self.assertEqual(self.root_logger.level, logging.DEBUG)
            self.assertEqual(len(self.root_logger.handlers), 2)
            self.assertTrue(any(isinstance(handler, logging.FileHandler) for handler in self.root_logger.handlers))
            self.assertTrue(any(isinstance(handler, logging.StreamHandler) and not isinstance(handler, logging.FileHandler) for handler in self.root_logger.handlers))

            log_dir = Path(temp_dir.name) / "logs"
            self.assertTrue(log_dir.exists())
            self.assertEqual(len(list(log_dir.glob("*.log"))), 1)
        finally:
            for handler in self.root_logger.handlers[:]:
                self.root_logger.removeHandler(handler)
                handler.close()
            temp_dir.cleanup()