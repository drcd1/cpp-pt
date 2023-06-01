using Plots

function importance(x)
    return x.*x.*cos.(x.*7).*cos.(x.*7).+0.1
end

function binary_search(cdf,ecs)

    ecs = cdf[end]*ecs
    if(ecs<=cdf[1])
        return 1
    end
    
    min = 2
    max = length(cdf)

    while max-min > 4

        mid = floor(Int,(max+min)/2)
        if ecs>cdf[mid]
            min = mid+1
        else
            max = mid
        end
        
    end

    for i = min:max
        if ecs<cdf[i]
            return i
        end
    end
    throw(ErrorException("Should not reach"))


end

function resample_binary(x,f)
    y=copy(x)
    cdf = copy(f)
    for i = 2:length(x)
        cdf[i]=cdf[i]+cdf[i-1]
    end

    for i = 1:length(y)
        ecs = rand()
        y[i] =x[binary_search(cdf,ecs)]
    end

    return y
end

function resample_two_steps(x,f)
    y = copy(x)
    #y = y[f.>rand(length(y))/length(y)*sum(f)]
    y = y[f.>rand(length(y))*maximum(f)]
    y2 = zeros(Float32, length(x)-length(y))
    for i=1:length(y2)
        y2[i] = y[floor(Int32,length(y)*rand())+1]
    end
    return [y;y2]
end

function histogram(x, intervals)
    a = zeros(length(intervals)-1)
    #assume intervals are regular
    iw = (intervals[end]-intervals[1])
    for i=1:length(x)
        helper = (x[i]-intervals[1])/iw*(length(intervals)-1);
        idx = floor(Int,helper) + 1
        a[idx] += 1
    end

    a = a./sum(a);
    return a;    

end

function integral(f,intervals)
    n=100
    a = zeros(length(intervals)-1)
    for i = 1:length(intervals)-1
        for i=1:n
            t = rand()
            a[i]+=f(intervals[i]*(1.0-t)+intervals[i+1]*t)/n
        end
    end

    a=a./sum(a)
    return a

end




#resamle binary converges faster?
#do they both converge?
#if new samples are resampled
#if new samples are 

function main()
    println("Hello world!")

    x = [0:9999;].*0.0001
    f = importance(x)
    y1 = resample_binary(x,f)
    y2 = resample_two_steps(x,f)

    intervals = [0:100;]*0.01;


    h1 = histogram(x,intervals)
    h2 = histogram(y1,intervals)
    h3 = histogram(y2,intervals)
    h4 = integral(importance,intervals)
    
    x_pts = (intervals[1:end-1]+intervals[2:end]).*0.5
    default(show=true)
    return scatter(x_pts,[h1,h2,h3,h4])
    


end



main()