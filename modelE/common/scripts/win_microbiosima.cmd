@echo off
setlocal
if ""=="%MICROBIOSIMA%" set MICROBIOSIMA=%~dp0%..
set MICROBIOSIMA_LIB=%MICROBIOSIMA%\lib
java -Xms64m -Xmx256m -Djava.library.path="%MICROBIOSIMA_LIB%" -jar "%MICROBIOSIMA_LIB%/microbiosima.jar" %*
