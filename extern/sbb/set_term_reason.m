function r = set_term_reason(t)

switch t
  case 1
    r = 'Exceeded time limit';
  case 2
    r = 'Relative change in x small enough';
  case 3
    r = 'Relative change in objvalue small enough';
  case 4
    r = 'Maximum number of iterations reached';
  case 5
    r = '|x_t+1 - x_t|=0 or |grad_t+1 - grad_t| < 1e-9';
  case 6
    r = 'Line search faild';
  case 7
    r = '|| proj(grad) ||_inf < opt.tolk';
 case 8
    r = '|| grad ||_inf < opt.tolg';
 case 10
    r = '|| proj(x - grad) - x ||_inf < opt.tolp';
  case 11
    r = 'norm(Ag) is too small';
  case 100
    r = 'The active set converged';
  otherwise
    r = 'Undefined';
end
