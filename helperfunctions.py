
def linspace(start, stop, n):
    if n == 1:
        return start
    else:
        output = []
        h = (stop - start) / (n - 1)
        for i in range(n):
            output.append(start + h * i)
        return output

