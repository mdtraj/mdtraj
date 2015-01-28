xcopy %RECIPE_DIR%\\..\\.. %SRC_DIR% /e /h /Y
"%PYTHON%" setup.py install
if errorlevel 1 exit 1
