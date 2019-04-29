import os
from urllib.parse import urlparse

from idseq_dag.util.s3 import check_s3_presence_for_file_list, fetch_from_s3, check_s3_presence, upload_with_retries
from shutil import copy


# Collection of classes and utility functions to work with files and directories including those that are found
# in S3 and possibly other key-value stores.

class Path():
    """
    Abstract path formed from url provided to the constructor. Unlike generic url we check that only specific schema
    are allowed. Currently we allow only s3:// and file://
    """
    def __init__(self, url_str):
        self.url = urlparse(url_str)
        if self.url.scheme == '':
            raise ValueError("Missing scheme")
        elif self.url.scheme != 's3' and self.url.scheme != 'file':
            raise ValueError("Unrecognized scheme")
        self.path = os.path.join("/" + self.url.netloc + "/", self.url.path[1:].rstrip('/'))

class File(Path):
    """
    Object that represents file.
    """
    def __init__(self, url_str):
        super(File, self).__init__(url_str)
        # We verify that object referred by url_string does not exist or is not a file only in case were scheme
        # is file:// i.e it refers to object in the actual filesystem. Later we can add support for s3:// or
        # any there key-value store that we will choose to support later.
        if self.url.scheme == 'file' and os.path.exists(self.path) \
                and not os.path.isfile(self.path):
            raise RuntimeError("Local " + self.path + " is a directory while expected to be file or non-existent.")

    def withNewPath(self, newPath):
        """
        Create new instance of File object with new path replacing old one.
        """
        return File(self.url.scheme + ':' + newPath)

    def copyTo(self, destination_dir):
        """
        Copy this file to destination directory in the local filesystem
        :return: Full path to the copied file in the local filesystem
        """
        os.makedirs(destination_dir, exist_ok=True)
        if self.url.scheme == "s3":
            return fetch_from_s3(str(self.url), destination_dir, allow_s3mi=True)
        elif self.url.scheme == "file":
            return str(copy(self.path, destination_dir))

    def copyFrom(self, source_file):
        """
        Copy source file into destination indicated by this File object
        :param source_file:
        :return:
        """
        if self.url.scheme == "s3":
            upload_with_retries(source_file, str(self.url.geturl()))
        elif self.url.scheme == "file":
            os.makedirs(os.path.dirname(self.path), exist_ok=True)
            copy(source_file, self.path)

    def exist(self):
        if self.url.scheme == "s3":
            return check_s3_presence(self.url.geturl())
        elif self.url.scheme == "file":
            return os.path.exists(self.path) and os.path.isfile(self.path)

class Dir(Path):
    """
    Object that represents directory.
    """
    def __init__(self, url_str):
        super(Dir, self).__init__(url_str)
        if self.url.scheme == 'file' and os.path.exists(self.path) and not os.path.isdir(self.path):
            raise RuntimeError("Local " + self.path + " is not a directory while expected to be dir or non-existent.")

    def withNewPath(self, newPath):
        """
        Create new instance of Dir object with new path replacing old one but the same scheme.
        """
        return Dir(self.url.scheme + ':' + newPath)

    def check_files_exist(self, file_list):
        """
        Verify that files provided in the file_list are present under this directory returning True if all exist
        and False otherwise.
        """
        if self.url.scheme == "s3":
            return check_s3_presence_for_file_list(self.url.geturl(), file_list)
        else:
            for file_name in file_list:
                full_file_name = os.path.join(self.path, file_name)
                if not os.path.exists(full_file_name) or not os.path.isfile(full_file_name):
                    return False
            return True
