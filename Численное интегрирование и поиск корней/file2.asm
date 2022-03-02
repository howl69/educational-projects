section .data
    f1_const1 dq 2.
    f2_const1 dq -2.
    f2_const2 dq 8.
    f3_const1 dq -5.
    df2_const1 dq -2.
    df3_const1 dq 5.
    
section .text
global f1
f1:
    push ebp
    mov ebp, esp
    finit
    fld qword [ebp + 8]
    fldl2e
    fmulp
    fld ST0
    fld1
    fxch ST0, ST2
    frndint
    fsub ST1, ST0
    fxch ST1, ST0
    f2xm1
    fadd ST0, ST2
    fscale
    fld qword[f1_const1]
    faddp
    mov esp, ebp
    pop ebp
    ret
    
global f2
f2:
    push ebp
    mov ebp, esp
    finit
    fld qword[ebp+8]
    fld qword[f2_const1]
    fmulp
    fld qword[f2_const2]
    faddp
    pop ebp
    ret
    
global f3
f3:
    push ebp
    mov ebp, esp
    finit
    fld qword[f3_const1]
    fld qword[ebp+8]
    fdivp
    pop ebp
    ret

global df1 ;exp(x)
df1:
    push ebp
    mov ebp, esp
    finit
    fld qword [ebp + 8]
    fldl2e
    fmulp
    fld ST0
    fld1
    fxch ST0, ST2
    frndint
    fsub ST1, ST0
    fxch ST1, ST0
    f2xm1
    fadd ST0, ST2
    fscale
    mov esp, ebp
    pop ebp
    ret

global df2; -2
df2:
    push ebp
    mov ebp, esp
    finit
    fld qword[df2_const1]
    pop ebp
    ret

global df3;5/(x*x)
df3:
    push ebp
    mov ebp, esp
    finit
    fld qword[df3_const1]
    fld qword[ebp+8]
    fdivp
    fld qword[ebp+8]
    fdivp
    pop ebp
    ret

global f4; sin(x)
f4:
    push ebp
    mov ebp, esp
    finit
    fld qword [ebp + 8]
    fsin    
    mov esp, ebp
    pop ebp
    ret
    
global f5; cos(x)
f5:
    push ebp
    mov ebp, esp
    finit
    fld qword [ebp + 8]
    fcos  
    mov esp, ebp
    pop ebp
    ret
    
global f6; 2^x + 2
f6:
    push ebp
    mov ebp, esp
    finit
    fld qword [ebp + 8]
    fld ST0
    fld1
    fxch ST0, ST2
    frndint
    fsub ST1, ST0
    fxch ST1, ST0
    f2xm1
    fadd ST0, ST2
    fscale
    fadd ST0, ST2 
    fadd ST0, ST2
    mov esp, ebp
    pop ebp
    ret
    
global f7; 4^x
f7:
    push ebp
    mov ebp, esp
    finit
    fld1
    fld qword [ebp + 8]
    fscale
    fld ST0
    frndint
    fsub ST1, ST0
    fxch ST1, ST0
    f2xm1
    fadd ST0, ST2
    fscale
    fsub ST0, ST2
    fld1
    fadd ST0, ST1
    mov esp, ebp
    pop ebp
    ret
    
global f8
f8: ;log2(x - 1)
    push ebp
    mov ebp, esp
    finit
    fld1
    fchs
    fld qword [ebp + 8]
    fadd ST0, ST1
    fyl2x
    fchs
    mov esp, ebp
    pop ebp
    ret
    
global f9
f9: ; -log2(1 + x)
    push ebp
    mov ebp, esp
    finit
    fld1
    fld qword [ebp + 8]
    fyl2xp1
    fchs
    mov esp, ebp
    pop ebp
    ret
    