##
# File:  PdbxReader.py
# Date:  2012-01-09  Jdw  Adapted from PdbxParser
#
# Updates:
#
# 2012-01-09 - (jdw) Separate reader and writer classes.
#
# 2012-09-02 - (jdw)  Revise tokenizer to better handle embedded quoting.
#
# 2025-05-15 - (cb) Changes docstrings to mdtraj format, exception handling.
#                   Replaces import directory from OpenMM.
#                   Removes unused __tokenizerOrg
##
"""
PDBx/mmCIF dictionary and data file parser.

Acknowledgements
---------------
The tokenizer used in this module is modeled after the clever parser design
used in the PyMMLIB package.

PyMMLib Development Group
Authors: Ethan Merritt (merritt@u.washington.ed) & Jay Painter (jay.painter@gmail.com)
See: http://pymmlib.sourceforge.net/
"""

import re

from .PdbxContainers import DataCategory, DataContainer, DefinitionContainer


class PdbxError(Exception):
    """
    Class for catch general errors
    """

    pass


class SyntaxError(Exception):
    """
    Exception raised for syntax errors in the mmCIF file.

    Parameters
    ----------
    lineNumber : int
        The line number where the error occurred.
    text : str
        Description of the syntax error.
    """

    def __init__(self, lineNumber, text):
        super().__init__()
        self.lineNumber = lineNumber
        self.text = text

    def __str__(self):
        return f"%ERROR - [at line: {self.lineNumber}] {self.text}"


class PdbxReader:
    """
    PDBx reader for data files and dictionaries.

    Parameters
    ----------
    ifh : file-like object
        Input file handle returned by open().
    """

    def __init__(self, ifh):
        """ifh - input file handle returned by open()"""
        #
        self.__curLineNumber = 0
        self.__ifh = ifh
        self.__stateDict = {
            "data": "ST_DATA_CONTAINER",
            "loop": "ST_TABLE",
            "global": "ST_GLOBAL_CONTAINER",
            "save": "ST_DEFINITION",
            "stop": "ST_STOP",
        }

    def read(self, containerList):
        """
        Parse the input file and append data and definition containers to containerList.

        Parameters
        ----------
        containerList : list
            List to which the parsed data and definition objects will be appended.

        Raises
        ------
        PdbxError
            If an unexpected error occurs during parsing.
        """
        self.__curLineNumber = 0
        try:
            self.__parser(self.__tokenizer(self.__ifh), containerList)
        except StopIteration:
            pass
        except RuntimeError as err:
            if "StopIteration" not in str(err):
                raise
        else:
            raise PdbxError()

    def __syntaxError(self, errText):
        raise SyntaxError(self.__curLineNumber, errText)

    def __getContainerName(self, inWord):
        """
        Returns the name of the data_ or save_ container

        Parameters
        ----------
        inWord : str
            A reserved word starting with 'data_' or 'save_'.

        Returns
        -------
        str
            The container name.
        """
        return str(inWord[5:]).strip()

    def __getState(self, inWord):
        """
        Identifies reserved syntax elements and assigns an associated state.

        Parameters
        ----------
        inWord : str
            Token representing a reserved syntax element.

        Returns
        -------
        tuple
            A tuple (reserved_word, state) where reserved_word is one of the CIF syntax elements
            (e.g., data_, loop_, global_, save_, stop_) and state is the corresponding parser state
            required to process the next section.
        """
        i = inWord.find("_")
        if i == -1:
            return None, "ST_UNKNOWN"

        try:
            rWord = inWord[:i].lower()
            return rWord, self.__stateDict[rWord]
        except KeyError:
            return None, "ST_UNKNOWN"

    def __parser(self, tokenizer, containerList):
        """
        Parse PDBx data files and dictionaries.

        This method iterates over the tokens from the tokenizer and populates the
        containerList with data and definition objects.

        Parameters
        ----------
        tokenizer : generator
            A reentrant method recognizing data item names (_category.attribute)
            quoted strings (single, double and multi-line semi-colon delimited),
            and unquoted strings.
        containerList : list
            Container for data and definition objects parsed from the input file.

        Raises
        ------
        SyntaxError
            If any syntax error is encountered during parsing.
        """
        # Working container - data or definition
        curContainer = None
        #
        # Working category container
        categoryIndex = {}
        curCategory = None
        #
        curRow = None
        state = None

        # Find the first reserved word and begin capturing data.
        #
        while True:
            curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
            if curWord is None:
                continue
            reservedWord, state = self.__getState(curWord)
            if reservedWord is not None:
                break

        while True:
            #
            #  Set the current state  -
            #
            #  At this point in the processing cycle we are expecting a token containing
            #  either a '_category.attribute'  or a reserved word.
            #
            if curCatName is not None:
                state = "ST_KEY_VALUE_PAIR"
            elif curWord is not None:
                reservedWord, state = self.__getState(curWord)
            else:
                self.__syntaxError("Miscellaneous syntax error")
                return

            #
            # Process  _category.attribute  value assignments
            #
            if state == "ST_KEY_VALUE_PAIR":
                try:
                    curCategory = categoryIndex[curCatName]
                except KeyError:
                    # A new category is encountered - create a container and add a row
                    curCategory = categoryIndex[curCatName] = DataCategory(curCatName)

                    try:
                        curContainer.append(curCategory)
                    except AttributeError:
                        self.__syntaxError("Category cannot be added to  data_ block")
                        return

                    curRow = []
                    curCategory.append(curRow)
                else:
                    # Recover the existing row from the category
                    try:
                        curRow = curCategory[0]
                    except IndexError:
                        self.__syntaxError(
                            "Internal index error accessing category data",
                        )
                        return

                # Check for duplicate attributes and add attribute to table.
                if curAttName in curCategory.getAttributeList():
                    self.__syntaxError("Duplicate attribute encountered in category")
                    return
                else:
                    curCategory.appendAttribute(curAttName)

                # Get the data for this attribute from the next token
                tCat, tAtt, curQuotedString, curWord = next(tokenizer)

                if tCat is not None or (curQuotedString is None and curWord is None):
                    self.__syntaxError(
                        f"Missing data for item _{curCatName}.{curAttName}",
                    )

                if curWord is not None:
                    #
                    # Validation check token for misplaced reserved words  -
                    #
                    reservedWord, state = self.__getState(curWord)
                    if reservedWord is not None:
                        self.__syntaxError(
                            f"Unexpected reserved word: {reservedWord}",
                        )

                    curRow.append(curWord)

                elif curQuotedString is not None:
                    curRow.append(curQuotedString)

                else:
                    self.__syntaxError("Missing value in item-value pair")

                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)
                continue

            #
            # Process a loop_ declaration and associated data -
            #
            elif state == "ST_TABLE":
                # The category name in the next curCatName,curAttName pair
                #    defines the name of the category container.
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

                if curCatName is None or curAttName is None:
                    self.__syntaxError("Unexpected token in loop_ declaration")
                    return

                # Check for a previous category declaration.
                if curCatName in categoryIndex:
                    self.__syntaxError("Duplicate category declaration in loop_")
                    return

                curCategory = DataCategory(curCatName)

                try:
                    curContainer.append(curCategory)
                except AttributeError:
                    self.__syntaxError(
                        "loop_ declaration outside of data_ block or save_ frame",
                    )
                    return

                curCategory.appendAttribute(curAttName)

                # Read the rest of the loop_ declaration
                while True:
                    curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

                    if curCatName is None:
                        break

                    if curCatName != curCategory.getName():
                        self.__syntaxError("Changed category name in loop_ declaration")
                        return

                    curCategory.appendAttribute(curAttName)

                # If the next token is a 'word', check it for any reserved words -
                if curWord is not None:
                    reservedWord, state = self.__getState(curWord)
                    if reservedWord is not None:
                        if reservedWord == "stop":
                            return
                        else:
                            self.__syntaxError(
                                f"Unexpected reserved word after loop declaration: {reservedWord}",
                            )

                # Read the table of data for this loop_ -
                while True:
                    curRow = []
                    curCategory.append(curRow)

                    for tAtt in curCategory.getAttributeList():
                        if curWord is not None:
                            curRow.append(curWord)
                        elif curQuotedString is not None:
                            curRow.append(curQuotedString)

                        curCatName, curAttName, curQuotedString, curWord = next(
                            tokenizer,
                        )

                    # loop_ data processing ends if -

                    # A new _category.attribute is encountered
                    if curCatName is not None:
                        break

                    # A reserved word is encountered
                    if curWord is not None:
                        reservedWord, state = self.__getState(curWord)
                        if reservedWord is not None:
                            break

                continue

            elif state == "ST_DEFINITION":
                # Ignore trailing unnamed saveframe delimiters e.g. 'save_'
                sName = self.__getContainerName(curWord)
                if len(sName) > 0:
                    curContainer = DefinitionContainer(sName)
                    containerList.append(curContainer)
                    categoryIndex = {}
                    curCategory = None

                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

            elif state == "ST_DATA_CONTAINER":
                #
                dName = self.__getContainerName(curWord)
                if len(dName) == 0:
                    dName = "unidentified"
                curContainer = DataContainer(dName)
                containerList.append(curContainer)
                categoryIndex = {}
                curCategory = None
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

            elif state == "ST_STOP":
                return
            elif state == "ST_GLOBAL":
                curContainer = DataContainer("blank-global")
                curContainer.setGlobal()
                containerList.append(curContainer)
                categoryIndex = {}
                curCategory = None
                curCatName, curAttName, curQuotedString, curWord = next(tokenizer)

            elif state == "ST_UNKNOWN":
                self.__syntaxError(f"Unrecogized syntax element: {str(curWord)}")
                return

    def __tokenizer(self, ifh):
        """
        Tokenize the mmCIF file.

        This generator method yields tokens from the input file. Each yielded tuple
        has the following structure:
            (category name, attribute name, quoted string, unquoted word)

        Parameters
        ----------
        ifh : file-like object
            The input file handle.

        Yields
        ------
        tuple
            A tuple containing:
            - category name : str or None
            - attribute name : str or None
            - quoted string : str or None
            - unquoted word : str or None
        """
        #
        # Regex definition for mmCIF syntax - semi-colon delimited strings are handled
        #                                     outside of this regex.
        mmcifRe = re.compile(
            r"(?:"
            r"(?:_(.+?)[.](\S+))"
            "|"  # _category.attribute
            r"(?:['](.*?)(?:[']\s|[']$))"
            "|"  # single quoted strings
            r"(?:[\"](.*?)(?:[\"]\s|[\"]$))"
            "|"  # double quoted strings
            r"(?:\s*#.*$)"
            "|"  # comments (dumped)
            r"(\S+)"  # unquoted words
            r")",
        )

        fileIter = iter(ifh)

        ## Tokenizer loop begins here ---
        while True:
            line = next(fileIter)
            self.__curLineNumber += 1

            # Dump comments
            if line.startswith("#"):
                continue

            # Gobble up the entire semi-colon/multi-line delimited string and
            #    and stuff this into the string slot in the return tuple
            #
            if line.startswith(";"):
                mlString = [line[1:]]
                while True:
                    line = next(fileIter)
                    self.__curLineNumber += 1
                    if line.startswith(";"):
                        break
                    mlString.append(line)

                # remove trailing new-line that is part of the \n; delimiter
                mlString[-1] = mlString[-1].rstrip()
                #
                yield (None, None, "".join(mlString), None)
                #
                # Need to process the remainder of the current line -
                line = line[1:]
                # continue

            # Apply regex to the current line consolidate the single/double
            # quoted within the quoted string category
            for it in mmcifRe.finditer(line):
                tgroups = it.groups()
                if tgroups != (None, None, None, None, None):
                    if tgroups[2] is not None:
                        qs = tgroups[2]
                    elif tgroups[3] is not None:
                        qs = tgroups[3]
                    else:
                        qs = None
                    groups = (tgroups[0], tgroups[1], qs, tgroups[4])
                    yield groups
