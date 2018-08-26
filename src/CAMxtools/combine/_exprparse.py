
"""
@author: 2017-01-27 jjung 

Most structure of this module is from,
Section 2.19 of "Python Cookbook 3rd edition" by David Beazley and Brian K. Jones in O'Reilly.

The original example is supposed to use a variable that has a single value such as 1.0 or 2.0. I develop the original module to pass arrays (Numpy arrays can be also used).

"""
import re
import collections

# token spec
NAME   = r'(?P<NAME>[a-zA-Z_][a-zA-Z_0-9\[\]]*)'
NUM    = r'(?P<NUM>[\d\.]+)'
PLUS   = r'(?P<PLUS>\+)'
MINUS  = r'(?P<MINUS>\-)'
TIMES  = r'(?P<TIMES>\*)'
DIVIDE = r'(?P<DIVIDE>/)'
LPAREN = r'(?P<LPAREN>\()'
RPAREN = r'(?P<RPAREN>\))'
WS     = r'(?P<WS>\s+)'

master_pat = re.compile('|'.join([NUM, PLUS, MINUS, TIMES,
                                  DIVIDE, LPAREN, RPAREN, WS, NAME]))
# token
Token = collections.namedtuple('Token',['type','value'])

def generate_tokens(text):
    scanner = master_pat.scanner(text)
    for m in iter(scanner.match, None):
        tok = Token(m.lastgroup, m.group())
        if tok.type !='WS':
            yield tok

# parser
class ExpressionEvaluator:
    '''
    recursive parser, all methods follow one grammar rule.
    Currently receive lookahead token and use ._accept() to test
    When input exactly matches and ignore next tokens,
    use ._expect() (If not matched, return SyntaxError.)
    '''

    def __init__(self,var):
        self.var = var
    def parse(self,text):
        self.tokens = generate_tokens(text)
        self.tok = None              # Use the last symbol
        self.nexttok = None          # Tokenize the next symbol
        self._advance()              # Call the first lookahead token
        return self.expr()

    def _advance(self):
        'Advance one token ahead'
        self.tok, self.nexttok = self.nexttok, next(self.tokens, None)

    def _accept(self,toktype):
        'Test and consume the next token if it matches toktype'
        if self.nexttok and self.nexttok.type == toktype:
            self._advance()
            return True
        else:
            return False

    def _expect(self,toktype):
        'Consume next token if it matches toktype or raise SyntaxError'
        if not self._accept(toktype):
            raise SyntaxError('Expect ' + toktype)

    # Grammar rule

    def expr(self):
        "expression ::= term { ('+'|'-') term }*"

        exprval = self.term()
        while self._accept('PLUS') or self._accept('MINUS'):
            op = self.tok.type
            right = self.term()
            if op == 'PLUS':
                exprval += right
            elif op == 'MINUS':
                exprval -= right
        return exprval

    def term(self):
        "term ::= factor { ('*'|'/') factor }*"

        termval = self.factor()
        while self._accept('TIMES') or self._accept('DIVIDE'):
            op = self.tok.type
            right = self.factor()
            if op == 'TIMES':
                termval *= right
            elif op == 'DIVIDE':
                termval /= right
        return termval

    def factor(self):
        "factor ::= { ('NUM'|'NAME') | ( expr )}"

        if self._accept('NUM'):
            return float(self.tok.value)
        elif self._accept('NAME'):
            return 1.0*self.var[self.tok.value]
        elif self._accept('LPAREN'):
            exprval = self.expr()
            self._expect('RPAREN')
            return exprval
        else:
            raise SyntaxError('Expected NUMBER or LPAREN')

class ExpressionTreeBuilder(ExpressionEvaluator):
    def expr(self):
        "expression ::= term { ('+'|'-') term }"

        exprval = self.term()
        while self._accept('PLUS') or self._accept('MINUS'):
            op = self.tok.type
            right = self.term()
            if op == 'PLUS':
                exprval = ('+', exprval, right)
            elif op == 'MINUS':
                exprval = ('-', exprval, right)
        return exprval

    def term(self):
        "term ::= factor { ('*'|'/') factor }"

        termval = self.factor()
        while self._accept('TIMES') or self._accept('DIVIDE'):
            op = self.tok.type
            right = self.factor()
            if op == 'TIMES':
                termval = ('*', termval, right)
            elif op == 'DIVIDE':
                termval = ('/', termval, right)
        return termval

    def factor(self):
        "factor ::= { ('NUM'|'NAME') | ( expr )}"

        if self._accept('NUM') or self._accept('NAME'):
            #return int(self.tok.value)
            return self.tok.value
        elif self._accept('LPAREN'):
            exprval = self.expr()
            self._expect('RPAREN')
            return exprval
        else:
            raise SyntaxError('Expected NUMBER or LPAREN')
