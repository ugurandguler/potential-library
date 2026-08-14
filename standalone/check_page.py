#!/usr/bin/env python3
"""
Does the built page's script still close every string it opens?

There is no JavaScript engine on this machine, and a page whose template
literals do not balance fails in the one way that is hardest to notice: it
does not error visibly, it renders an element panel as empty or truncated.
Everything downstream - the tables, the plots, the notes - is generated at
runtime from those templates, so a single unclosed backtick silently deletes
a section that the build reported as written.

This is not a parser.  It tracks the four things that can swallow the rest of
a file and nothing else: backtick templates, the ${ } holes inside them,
ordinary quoted strings, and comments.

Comments are the part that has to be right rather than approximately right.
The source is full of prose - "ruthenium's page", "MP's unrounded tensor",
"doesn't" - and an apostrophe inside a comment opens a string that runs to the
next apostrophe hundreds of lines away.  A first version of this check skipped
comments, reported the page broken, and sent a long search after a fault that
was in the checker.  It is kept honest here by being run against a page known
to render before being trusted about one that has changed.

    python check_page.py                    # the current page
    python check_page.py a.html b.html      # compare, e.g. against a backup
"""
import io
import re
import sys

DIV = set("(,=:[!&|?{};+-*%~^<>")      # after these, / starts a regex


def scan(path):
    """(unclosed, description) - unclosed is 0 when the script balances"""
    html = io.open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"<script>([\s\S]*)</script>", html)
    if not m:
        return -1, "script not found"
    js = m.group(1)
    #  Each ${ } hole carries its OWN brace counter.  Sharing one counter
    #  across nested holes makes an inner hole's closing brace cancel an outer
    #  hole's opening one, and the stack then climbs and never comes back
    #  down.  That is not a subtle wrongness: it reported a page as broken
    #  from its first nested hole onwards, which is the second line of any
    #  template here, and the search for the imaginary fault went as far as
    #  removing a whole section before the checker itself was suspected.
    #  A frame is ("tpl", 0) or ("hole", braces seen so far).
    stack, i, n, q = [], 0, len(js), None
    prev = ""                       # last significant character
    while i < n:
        ch = js[i]
        intpl = bool(stack) and stack[-1][0] == "tpl"
        if q:
            if ch == "\\":
                i += 2
                continue
            if ch == q:
                q = None
            i += 1
            continue
        #  comments, but only where code can appear - inside a template body
        #  a slash is text
        if not intpl and ch == "/" and i + 1 < n:
            if js[i + 1] == "/":
                i = js.find("\n", i)
                i = n if i < 0 else i + 1
                continue
            if js[i + 1] == "*":
                j = js.find("*/", i + 2)
                i = n if j < 0 else j + 2
                continue
            #  a regex literal, which may contain quotes and backticks
            if prev in DIV or prev == "":
                j, esc, cls = i + 1, False, False
                while j < n:
                    c = js[j]
                    if esc:
                        esc = False
                    elif c == "\\":
                        esc = True
                    elif c == "[":
                        cls = True
                    elif c == "]":
                        cls = False
                    elif c == "/" and not cls:
                        break
                    elif c == "\n":
                        j = -1
                        break
                    j += 1
                if j > 0:
                    i = j + 1
                    prev = "/"
                    continue
        if ch == "\\":
            i += 2
            continue
        if ch in "'\"":
            if intpl:
                i += 1
                continue
            q = ch
            i += 1
            continue
        if ch == "`":
            if intpl:
                stack.pop()
            else:
                stack.append(["tpl", 0])
            i += 1
            prev = "`"
            continue
        if intpl:
            if ch == "$" and i + 1 < n and js[i + 1] == "{":
                stack.append(["hole", 0])
                i += 2
            else:
                i += 1
            continue
        if stack and stack[-1][0] == "hole":
            if ch == "{":
                stack[-1][1] += 1
            elif ch == "}":
                if stack[-1][1] == 0:
                    stack.pop()
                else:
                    stack[-1][1] -= 1
        if not ch.isspace():
            prev = ch
        i += 1
    if q:
        return len(stack) + 1, f"kapanmamis dizgi ({q})"
    return len(stack), (f"kapanmamis: {[f[0] for f in stack]}"
                        if stack else "dengeli")


def main():
    paths = sys.argv[1:] or ["potential.html"]
    bad = 0
    for p in paths:
        n, why = scan(p)
        print(f"{p[:44]:46s} {'TAMAM' if n == 0 else 'BOZUK'}  {why[:60]}")
        bad += (n != 0)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
