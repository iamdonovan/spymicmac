spymicmac.ee_tools
=================================

.. note::
    Before using spymicmac.ee_tools, you will need to run ee.Authenticate(), in order to authorize the access needed
    by EarthEngine. To do this from the python terminal, run the following:
    ::

        import ee
        ee.Authenticate()

    This will open a browser window and prompt you to select your Google account for authentication. It will generate
    an authorization code which you should then copy and paste into the prompt in the python terminal.

    After you have successfully authenticated, you should be able to use spymicmac.ee_tools normally (note that you
    will need to explicitly import ``ee_tools`` - it is not loaded by ``spymicmac.__init__.py``, in order to prevent
    authorization errors):
    ::

        from spymicmac import ee_tools

    You won't need to do this every time you use :py:meth:`spymicmac.ee_tools`, but you may need to re-authorize from time
    to time.

.. automodule:: spymicmac.ee_tools
    :members: