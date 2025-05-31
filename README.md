# Vendor Branch: port/openmm-pdbx

This branch serves as a **vendor branch** for importing and tracking a specific set of files from the [OpenMM](https://github.com/openmm/openmm) repository into [MDTraj](https://github.com/mdtraj/mdtraj).
The main function of this branch is to allow changes to be tracked as parallel commits are done in the MDTraj codebase and in the OpenMM codebase.
It was created as an **orphan branch** to isolate the vendor files from the rest of the MDTraj codebase:

```bash
git checkout --orphan port/openmm-pdbx
````

### Tracked Files from OpenMM

The following files are tracked from OpenMM (typically from the `master` branch):

* `openmm/wrappers/python/openmm/app/internal/pdbx/*`
* `openmm/wrappers/python/openmm/app/pdbxfile.py`
* `openmm/wrappers/python/tests/TestPdbxFile.py`

These files were added directly to the initial commit of this branch and are preserved as-is, without any MDTraj-specific changes.

---

## How to Update This Branch

To update the vendor files when OpenMM makes changes:

1. **Check out this branch:**

   ```bash
   git checkout port/openmm-pdbx
   ```

2. **Remove the old vendor files:**

   ```bash
   git rm -r openmm/wrappers/python/openmm/app/internal/pdbx
   git rm    openmm/wrappers/python/openmm/app/pdbxfile.py
   git rm    openmm/wrappers/python/tests/TestPdbxFile.py
   ```

3. **Copy the new versions of the vendor files:**
   Replace the deleted files with their updated versions from the OpenMM repo (for example, by checking out the latest OpenMM and copying files over):

   ```bash
   cp -r ../openmm/wrappers/python/openmm/app/internal/pdbx .
   cp ../openmm/wrappers/python/openmm/app/pdbxfile.py .
   cp ../openmm/wrappers/python/tests/TestPdbxFile.py .
   ```

4. **Stage and commit the updates:**

   ```bash
   git add openmm/wrappers/python/openmm/app/internal/pdbx
   git add openmm/wrappers/python/openmm/app/pdbxfile.py
   git add openmm/wrappers/python/tests/TestPdbxFile.py
   git commit -m "Update vendor files from OpenMM commit <sha or tag>"
   ```

---

## Merging into MDTraj

Once the vendor branch is updated:

1. **Check out your integration branch** (e.g., `main`, `develop`, or a feature branch in MDTraj):

   ```bash
   git checkout main
   ```

2. **Merge or rebase the vendor branch:**

   ```bash
   git merge port/openmm-pdbx
   # or
   # git rebase port/openmm-pdbx
   ```

3. **Apply any further MDTraj-specific adaptations** in follow-up commits.

---

## Notes

* This branch is **not meant for development**. Do not modify vendor files directly here.
* All adaptations or fixes should be done in the integration branch after the vendor merge.

---

## Reference

OpenMM repository: [https://github.com/openmm/openmm](https://github.com/openmm/openmm)

MDTraj repository: [https://github.com/mdtraj/mdtraj](https://github.com/mdtraj/mdtraj)

## License

The OpenMM wrappers are distributed under the MIT License:

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject
to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

OpenMM uses the PDBx/mmCIF parser written by John Westbrook.  Which is distributed
under the Creative Commons Attribution 3.0 Unported license.  For details, see
https://creativecommons.org/licenses/by/3.0.  The library in OpenMM was modified 
for the openmm.app.internal module.
