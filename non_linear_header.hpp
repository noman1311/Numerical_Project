#include<bits/stdc++.h>
using namespace std;

double eqn(double x, vector<double>& coeff) 
{ 
    double sum{0.0};
    for(int i = 0; i<coeff.size(); i++){
        sum += (coeff[i] * pow(x,i));
    }
    return sum;
}

// Bisection method to find the root
void bisection() 
{
    int degree; cout<<"enter degree: "; cin>>degree;
    vector<double>coeff(degree+1);
    cout<<"enter the coefficients:\n";
    for(int i =degree; i>=0; i--)
        cin>>coeff[i];
    double a, b;
    double max_tolerance;
    cout<<"enter max tolerance:\n"; cin>>max_tolerance;
    cout<<"\nenter 2 brackets: "; cin>>a>>b;
    if(eqn(a, coeff) * eqn(b, coeff) > 0){
        cout<<"invalid brackets. :-(\nf(a)*f(b) must be less than 0.";
        return ;
    }

    int iteration = 0;  // Counter to track iterations
    while (true) 
    { 
        iteration++; 
        double x = (a + b) / 2;  // Midpoint of the interval
        cout << "Approximate Root = " << x << endl;

        // Check if we have found the exact root
        if (eqn(x,coeff) == 0.0) 
        { 
            cout << "Exact Root = " << x << endl; 
            cout << "Iterations = " << iteration << endl; 
            return; 
        }

        // Adjust the interval
        if (eqn(x,coeff) * eqn(a,coeff) < 0) 
            b = x; 
        else 
            a = x;  // No need for eqn(x) * eqn(b) < 0 check; it's implied

        // Check for convergence
        if (fabs(b - a) <= max_tolerance) 
        { 
            cout << "Converged Root = " << x << endl; 
            cout << "Iterations = " << iteration << endl; 
            return; 
        }
    } 
} 
void false_position() 
{
    int degree; cout<<"enter degree: "; cin>>degree;
    vector<double>coeff(degree+1);
    cout<<"enter the coefficients:\n";
    for(int i =degree; i>=0; i--)
        cin>>coeff[i];
    double a, b;
    double max_tolerance;
    cout<<"enter max tolerance:\n"; cin>>max_tolerance;
    cout<<"\nenter 2 brackets: "; cin>>a>>b;
    if(eqn(a, coeff) * eqn(b, coeff) > 0){
        cout<<"invalid brackets. :-(\nf(a)*f(b) must be less than 0.";
        return ;
    } 


    int iteration=0; 
    double old_x=a; 
    while(true) 
    { 
        iteration++; 
        double x=(a*eqn(b,coeff)-b*eqn(a,coeff))/(eqn(b,coeff)-eqn(a,coeff)); 
        cout<<"Approximate Root = "<<x<<endl; 
        if(eqn(x,coeff)==0.0) 
        { 
            cout<<"Exact Root = "<<x<<endl; 
            cout<<"Iteration = "<<iteration<<endl; 
            return; 
        } 
        else if(eqn(x,coeff)*eqn(a,coeff)<0) 
            b=x; 
        else if(eqn(x,coeff)*eqn(b,coeff)<0) 
            a=x; 
        if(fabs(x-old_x)<=max_tolerance) 
        { 
            cout<<"Converged Root = "<<x<<endl; 
            cout<<"Iteration = "<<iteration<<endl; 
            return; 
        } 
        old_x=x; 
    } 
}

double f(double x, vector<double>& coeff) 
{ 
    double sum{0.0};
    for(int i = 0; i<coeff.size(); i++){
        sum += (coeff[i] * pow(x,i));
    }
    return sum;
}

void secant()
{
    int degree; cout<<"enter degree: "; cin>>degree;
    vector<double>coeff(degree+1);
    cout<<"enter the coefficients:\n";
    for(int i =degree; i>=0; i--)
        cin>>coeff[i];
    double add, bdd;
    double max_tolerance; int max_itr;
    cout<<"enter max tolerance and iteration:\n"; cin>>max_tolerance>>max_itr;
    cout<<"\nenter 2 assumptions: ";
    double x1, x2; cin>>x1>>x2;
    if(abs(f(x1,coeff) - f(x2,coeff)) < max_tolerance){
        cout << "\n Wrong assumption !" << endl;
        return ;
    }

    bool flag = true;
    for(int i = 0; i<max_itr; i++){
        double x3 = (x1*f(x2,coeff) - x2*f(x1,coeff)) / (f(x2,coeff) - f(x1,coeff));
        double error = fabs(x3 - x2);
        if(f(x3,coeff) == 0.0){
            cout<<"exact root is: "<<x3<<" in "<<i+1<<"th iteration"<<endl;
            flag = false; break;
        }
        else if(error <= max_tolerance){
            cout<<"converged root is: "<<x3<<" in "<<i+1<<"th iteration"<<endl;
            flag = false; break;
        }
        else{
            cout<<"approximated root is: "<<x3<<" in "<<i+1<<"th iteration"<<endl;
            x1 = x2; x2 = x3;
        }
    }
    if(flag)
        cout<<"no root found in any iteration. :-( \ntry with more iterations. 0_0";
    return;
}

double df(double x, vector<double>& coeff)
{
    double sum{0.0}, ab = 1;
    for(int i = 1; i<coeff.size(); ++i){
        sum += ab*coeff[i]*i;
        ab *= x;
    }
    return sum;
}
void newton_raphson()
{
    int degree; cout<<"enter degree: "; cin>>degree;
    vector<double>coeff(degree+1);
    cout<<"enter the coefficients:\n";
    for(int i =degree; i>=0; i--)
        cin>>coeff[i];
    double a, b;
    double max_tolerance; int max_itr;
    cout<<"enter max tolerance and iteration:\n"; cin>>max_tolerance>>max_itr;
    cout<<"\nenter 1 assumptions: ";
    double x0; cin>>x0;

    bool flag = true;
    for(int i = 0; i<max_itr; i++){
        double x1 = x0 - (f(x0, coeff) / df(x0, coeff));
        double error = fabs(x1 - x0);
        if(f(x1, coeff) == 0.0){
            cout<<"exact root is: "<<x1<<" in "<<i+1<<"th iteration"<<endl;
            flag = false; break;
        }
        else if(error <= max_tolerance){
            cout<<"converged root is: "<<x1<<" in "<<i+1<<"th iteration"<<endl;
            flag = false; break;
        }
        else{
            cout<<"approximated root is: "<<x1<<" in "<<i+1<<"th iteration"<<endl;
            x0 = x1;
        }
    }
    if(flag)
        cout<<"no root found in any iteration. :-( \ntry with more iteraion.";
    return;
}
