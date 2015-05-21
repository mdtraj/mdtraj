xcopy %RECIPE_DIR%\\..\\.. %SRC_DIR% /e /h /Y /Q
"%PYTHON%" setup.py install
if errorlevel 1 exit 1
