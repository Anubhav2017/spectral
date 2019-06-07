#ifndef PRIORITYQUEUE_H
#define PRIORITYQUEUE_H

#include <QList>
#include <QPair>

template <typename T>
class PriorityQueue
{
public:
    PriorityQueue()
    {
    }

    void enqueue(T element, double value)
    {
        QPair<T, double> newEntry(element, value);
        if (m_list.isEmpty() || value >= m_list.last().second) {
            m_list.append(newEntry);
        } else {
            int index = 0;
            const int listSize = m_list.size();
            while (index < listSize && m_list[index].second < value)
                index++;
            m_list.insert(index, newEntry);
        }
    }

    T dequeue()
    {
        return m_list.takeFirst().first;
    }

    void refreshDistance(T element, double newDistance)
    {
        int index = 0;
        while (m_list[index].first != element)
            index++;

        m_list.removeAt(index);

        this->enqueue(element, newDistance);
    }

    bool isEmpty()
    {
        return m_list.isEmpty();
    }

    int getSize()
    {
        return m_list.size();
    }

    bool verifyQueue()
    {
        bool allGood = true;
        const int listSize = m_list.size();
        for (int i = 0; i < listSize - 1; i++)
            allGood = allGood && (m_list[i].second <= m_list[i+1].second);

        return allGood;
    }

    int find(T t) {
        for (int i = 0; i < m_list.size(); i++)
            if (m_list[i].first == t)
                return i;

        return -1;
    }

    int getValue(T t) {
        int index = this->find(t);
        if (index != -1)
            return m_list[this->find(t)].second;
        else
            return 0;
    }

private:
    QList<QPair<T, double> > m_list;
};

#endif // PRIORITYQUEUE_H
