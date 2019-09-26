m=Array{Float64}(undef,10)
x=Array{Float64}(undef,10)
p=Array{Float64}(undef,10)
t=Array{Float64}(undef,10)
h=Array{Float64}{undef,10}

function model!(F,x0)
    m=x0[1:10];x=x0[11:20];p=x0[21:30];t=x0[31:40];h=x0[41:50]
    x[1]=55;x[4]=60;m[1]=1;
    F[1]=m[1]-m[10]-m[6]
    F[2]=m[4]+m[7]-m[3]
    F[3]=m[8]-m[7]
    F[4]=m[9]-m[8]
    F[5]=m[10]-m[9]
    F[6]=m[2]-m[1]
    F[7]=m[6]-m[5]
    F[8]=m[5]-m[4]
    F[9]=m[3]-m[2]

    F[10]=m[1]*x[1]-m[6]*x[6]
    F[11]=m[4]*x[4]-m[3]*x[3]
    F[12]=x[8]-x[7]
    F[13]=x[9]-x[8]
    F[14]=x[10]-x[8]
    F[15]=x[2]-x[1]
    F[16]=x[6]-x[5]
    F[17]=x[5]-x[4]
    F[18]=x[3]-x[2]
    
end
